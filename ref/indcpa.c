#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "indcpa.h"
#include "poly.h"
#include "polyvec.h"
#include "rng.h"
#include "ntt.h"
#include "symmetric.h"

/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk
*              and the public seed used to generate the matrix A.
*
* Arguments:   uint8_t *r:          pointer to the output serialized public key
*              polyvec *pk:         pointer to the input public-key polyvec
*              const uint8_t *seed: pointer to the input public seed
**************************************************/
static void pack_pk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES])
{
  size_t i;
  polyvec_tobytes(r, pk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    r[i+KYBER_POLYVECBYTES] = seed[i];
}

/*************************************************
* Name:        unpack_pk
*
* Description: De-serialize public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk:             pointer to output public-key
*                                         polynomial vector
*              - uint8_t *seed:           pointer to output seed to generate
*                                         matrix A
*              - const uint8_t *packedpk: pointer to input serialized public key
**************************************************/
static void unpack_pk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES])
{
  size_t i;
  polyvec_frombytes(pk, packedpk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    seed[i] = packedpk[i+KYBER_POLYVECBYTES];
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key
*
* Arguments:   - uint8_t *r:  pointer to output serialized secret key
*              - polyvec *sk: pointer to input vector of polynomials (secret key)
**************************************************/
static void pack_sk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  polyvec_tobytes(r, sk);
}

/*************************************************
* Name:        unpack_sk
*
* Description: De-serialize the secret key;
*              inverse of pack_sk
*
* Arguments:   - polyvec *sk:             pointer to output vector of
*                                         polynomials (secret key)
*              - const uint8_t *packedsk: pointer to input serialized secret key
**************************************************/
static void unpack_sk(polyvec *sk,
                      const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec_frombytes(sk, packedsk);
}

/*************************************************
* Name:        pack_ciphertext
*
* Description: Serialize the ciphertext as concatenation of the
*              compressed and serialized vector of polynomials b
*              and the compressed and serialized polynomial v
*
* Arguments:   uint8_t *r: pointer to the output serialized ciphertext
*              poly *pk:   pointer to the input vector of polynomials b
*              poly *v:    pointer to the input polynomial v
**************************************************/
static void pack_ciphertext(uint8_t r[KYBER_INDCPA_BYTES],
                            polyvec *b,
                            poly *v)
{
  polyvec_compress(r, b);
  poly_compress(r+KYBER_POLYVECCOMPRESSEDBYTES, v);
}

/*************************************************
* Name:        unpack_ciphertext
*
* Description: De-serialize and decompress ciphertext from a byte array;
*              approximate inverse of pack_ciphertext
*
* Arguments:   - polyvec *b:       pointer to the output vector of polynomials b
*              - poly *v:          pointer to the output polynomial v
*              - const uint8_t *c: pointer to the input serialized ciphertext
**************************************************/
static void unpack_ciphertext(polyvec *b,
                              poly *v,
                              const uint8_t c[KYBER_INDCPA_BYTES])
{
  polyvec_decompress(b, c);
  poly_decompress(v, c+KYBER_POLYVECCOMPRESSEDBYTES);
}

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r:          pointer to output buffer
*              - unsigned int len:    requested number of 16-bit integers
*                                     (uniform mod q)
*              - const uint8_t *buf:  pointer to input buffer
*                                     (assumed to be uniform random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint16_t val0, val1;

  ctr = pos = 0;
  while(ctr < len && pos + 3 <= buflen) {
    val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
    val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;

    if(val0 < KYBER_Q)
      r[ctr++] = val0;
    if(ctr < len && val1 < KYBER_Q)
      r[ctr++] = val1;
  }

  return ctr;
}

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

/*************************************************
* Name:        gen_matrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *a:          pointer to ouptput matrix A
*              - const uint8_t *seed: pointer to input seed
*              - int transposed:      boolean deciding whether A or A^T
*                                     is generated
**************************************************/
#define GEN_MATRIX_NBLOCKS ((12*KYBER_N/8*(1 << 12)/KYBER_Q \
                             + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)
// Not static for benchmarking
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed)
{
  unsigned int ctr, i, j, k;
  unsigned int buflen, off;
  uint8_t buf[GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES+2];
  xof_state state;

  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      if(transposed)
        xof_absorb(&state, seed, i, j);
      else
        xof_absorb(&state, seed, j, i);

      xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
      buflen = GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES;
      ctr = rej_uniform(a[i].vec[j].coeffs, KYBER_N, buf, buflen);

      while(ctr < KYBER_N) {
        off = buflen % 3;
        for(k = 0; k < off; k++)
          buf[k] = buf[buflen - off + k];
        xof_squeezeblocks(buf + off, 1, &state);
        buflen = off + XOF_BLOCKBYTES;
        ctr += rej_uniform(a[i].vec[j].coeffs + ctr, KYBER_N - ctr, buf, buflen);
      }
    }
  }
}

#include "printIV.h"
//#define PRINT_KEYGEN_IV
//#define PRINT_ENC_IV
//#define PRINT_DEC_IV
#define PRINT_SEEDS
// #define PRINT_KEM_IV_CHECK

/*************************************************
* Name:        indcpa_keypair
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
                              (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
**************************************************/
void indcpa_keypair(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  unsigned int i;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  uint8_t nonce = 0;
  polyvec a[KYBER_K], e, pkpv, skpv;

  randombytes(buf, KYBER_SYMBYTES);
#ifdef PRINT_KEM_IV_CHECK
print_bytes("buf", buf, KYBER_SYMBYTES);
#endif
#ifdef PRINT_SEEDS
printf("------------------------------------------\n");
print_bytes("buf", buf, KYBER_SYMBYTES);
#endif


#ifdef PRINT_KEYGEN_IV
printf("------------------------------------------\n");
print_bytes("buf", buf, KYBER_SYMBYTES);
#endif
  hash_g(buf, buf, KYBER_SYMBYTES);
#ifdef PRINT_KEYGEN_IV
print_bytes("publicseed", publicseed, KYBER_SYMBYTES);
print_bytes("noiseseed", noiseseed, KYBER_SYMBYTES);
#endif

  gen_a(a, publicseed);
#ifdef PRINT_KEYGEN_IV
print_A(a);
#endif

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(&skpv.vec[i], noiseseed, nonce++);
#ifdef PRINT_KEYGEN_IV
print_polyvec("skpv",skpv);
#endif
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(&e.vec[i], noiseseed, nonce++);
#ifdef PRINT_KEYGEN_IV
print_polyvec("e",e);
#endif

  polyvec_ntt(&skpv);
#ifdef PRINT_KEYGEN_IV
print_polyvec("skpv_hat",skpv);
#endif
  polyvec_ntt(&e);
#ifdef PRINT_KEYGEN_IV
print_polyvec("e_hat",e);
#endif

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
    polyvec_pointwise_acc_montgomery(&pkpv.vec[i], &a[i], &skpv);
    poly_tomont(&pkpv.vec[i]);
  }
#ifdef PRINT_KEYGEN_IV
print_polyvec("A_hat@s_hat",pkpv);
#endif

  polyvec_add(&pkpv, &pkpv, &e);
  polyvec_reduce(&pkpv);
#ifdef PRINT_KEYGEN_IV
print_polyvec("t_hat",pkpv);
#endif

  pack_sk(sk, &skpv);
  pack_pk(pk, &pkpv, publicseed);
#ifdef PRINT_KEYGEN_IV
print_bytes("sk", sk, KYBER_INDCPA_SECRETKEYBYTES);
print_bytes("pk", pk, KYBER_INDCPA_PUBLICKEYBYTES);
printf("------------------------------------------\n");
#endif
}

/*************************************************
* Name:        indcpa_enc
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *c:           pointer to output ciphertext
*                                      (of length KYBER_INDCPA_BYTES bytes)
*              - const uint8_t *m:     pointer to input message
*                                      (of length KYBER_INDCPA_MSGBYTES bytes)
*              - const uint8_t *pk:    pointer to input public key
*                                      (of length KYBER_INDCPA_PUBLICKEYBYTES)
*              - const uint8_t *coins: pointer to input random coins
*                                      used as seed (of length KYBER_SYMBYTES)
*                                      to deterministically generate all
*                                      randomness
**************************************************/
void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES];
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, at[KYBER_K], bp;
  poly v, k, epp;

#ifdef PRINT_SEEDS
print_bytes("m"    , m    , KYBER_INDCPA_MSGBYTES);
print_bytes("coins", coins, KYBER_SYMBYTES);
#endif

#ifdef PRINT_ENC_IV
printf("------------------------------------------\n");
print_bytes("pk"   , pk   , KYBER_INDCPA_PUBLICKEYBYTES);
print_bytes("m"    , m    , KYBER_INDCPA_MSGBYTES);
print_bytes("coins", coins, KYBER_SYMBYTES);
#endif
  unpack_pk(&pkpv, seed, pk);
  poly_frommsg(&k, m);
  gen_at(at, seed);
#ifdef PRINT_ENC_IV
print_polyvec("t_hat" , pkpv);
print_bytes(  "rho"   , seed , KYBER_SYMBYTES);
print_poly(   "k"     , k   );
print_A(at);
#endif

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(sp.vec+i, coins, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta2(ep.vec+i, coins, nonce++);
  poly_getnoise_eta2(&epp, coins, nonce++);
#ifdef PRINT_ENC_IV
print_polyvec("sp" , sp );
print_polyvec("ep" , ep );
print_poly(   "epp", epp);
#endif

  polyvec_ntt(&sp);
#ifdef PRINT_ENC_IV
print_polyvec("sp_hat" , sp );
#endif
  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++)
    polyvec_pointwise_acc_montgomery(&bp.vec[i], &at[i], &sp);

  polyvec_pointwise_acc_montgomery(&v, &pkpv, &sp);
#ifdef PRINT_ENC_IV
print_polyvec("at_hat@sp_hat" , bp );
print_poly("t_hat@sp_hat" , v  );
#endif

  polyvec_invntt_tomont(&bp);
  poly_invntt_tomont(&v);
#ifdef PRINT_ENC_IV
print_polyvec("at@sp" , bp );
print_poly("t@sp" , v  );
#endif

  polyvec_add(&bp, &bp, &ep);
  poly_add(&v, &v, &epp);
#ifdef PRINT_ENC_IV
print_polyvec("at@sp+e1" , bp );
print_poly("t@sp+e2" , v  );
#endif
  poly_add(&v, &v, &k);
#ifdef PRINT_ENC_IV
print_poly("t@sp+e2+k" , v  );
#endif
  polyvec_reduce(&bp);
  poly_reduce(&v);
#ifdef PRINT_ENC_IV
print_polyvec("u" , bp );
print_poly("v" , v  );
#endif

  pack_ciphertext(c, &bp, &v);
#ifdef PRINT_ENC_IV
print_bytes("ciphertext" , c, KYBER_INDCPA_BYTES);
#endif
}

/*************************************************
* Name:        indcpa_dec
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *m:        pointer to output decrypted message
*                                   (of length KYBER_INDCPA_MSGBYTES)
*              - const uint8_t *c:  pointer to input ciphertext
*                                   (of length KYBER_INDCPA_BYTES)
*              - const uint8_t *sk: pointer to input secret key
*                                   (of length KYBER_INDCPA_SECRETKEYBYTES)
**************************************************/
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec bp, skpv;
  poly v, mp;

#ifdef PRINT_DEC_IV
printf("------------------------------------------\n");
print_bytes("c" , c, KYBER_INDCPA_BYTES);
print_bytes("sk" , sk, KYBER_INDCPA_SECRETKEYBYTES);
#endif
  unpack_ciphertext(&bp, &v, c);
  unpack_sk(&skpv, sk);
#ifdef PRINT_DEC_IV
print_polyvec("u" , bp );
print_poly("v" , v  );
print_polyvec("sk" , skpv );
#endif

  polyvec_ntt(&bp);
  polyvec_pointwise_acc_montgomery(&mp, &skpv, &bp);
#ifdef PRINT_DEC_IV
print_polyvec("u_hat" , bp );
print_poly("s_hat_u_hat" , mp);
#endif
  poly_invntt_tomont(&mp);
#ifdef PRINT_DEC_IV
print_poly("s_u" , mp);
#endif

  poly_sub(&mp, &v, &mp);
#ifdef PRINT_DEC_IV
print_poly("v_s_u" , mp);
#endif
  poly_reduce(&mp);
#ifdef PRINT_DEC_IV
print_poly("k" , mp);
#endif

  poly_tomsg(m, &mp);
#ifdef PRINT_DEC_IV
print_bytes("m" , m, KYBER_INDCPA_MSGBYTES);
#endif
}
