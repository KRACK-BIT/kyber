#ifndef _PRINTIV_H
#define _PRINTIV_H

#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "symmetric.h"

static void print_bytes(const char* label, const uint8_t* array, int len){
	printf("%s:", label);
	for(int i = 0; i < len; ++i){
		if(i%32 == 0){
			printf("\n");
		}
		printf("%02X ", array[i]);
	}
	printf("\n");
}

static void print_poly(const char* label, const poly p){
	printf("%s:", label);
	for(int i = 0; i < KYBER_N; ++i){
		if(i%32 == 0){
			printf("\n");
		}
		printf("%04d ", p.coeffs[i]);
	}
	printf("\n");
}

static void print_poly_hash(const char* label, const poly p){
	unsigned char buf[32];
	hash_h(buf, (char*)p.coeffs, KYBER_N*2);
	print_bytes(label, buf, 32);
}

static void print_polyvec(const char* label, const polyvec pv){
	char str[20];
	printf("%s:\n", label);
	for(int i = 0; i < KYBER_K; ++i){
		sprintf(str, "%s[%d]", label, i);
		print_poly(str, pv.vec[i]);
	}
	printf("\n");

	unsigned char buf[KYBER_N*KYBER_K*2];
	printf("H(%s):\n",label);
	for(int i = 0; i < KYBER_K; ++i){
		sprintf(str, "%s[%d]", label, i);
		print_poly_hash(str, pv.vec[i]);
	}
	printf("\n");

}

static void print_A(const polyvec* pv){
	char str[20];
	printf("A:\n");
	for(int i = 0; i < KYBER_K; ++i){
		for(int j = 0; j < KYBER_K; ++j){
			sprintf(str, "A[%d][%d]", i, j);
			print_poly(str, pv[i].vec[j]);
		}
	}
	printf("\n");
}

//int main()
//{
//	uint8_t buf[2*32] = {0};
//
//	for(int i = 0 ;i < 32; ++i){
//		buf[i] = (uint8_t)i;
//	}
//
//	print_bytes("test1", buf, 128);
//}

#endif
