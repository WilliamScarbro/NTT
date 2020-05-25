#ifndef NNT_H
#define NTT_H

void ntt(long* vec, long* dest, long N, long P, long root);
void inv_ntt(long* vec, long* dest, long N, long P, long root);
void printPoly(long* vec, long N);
//void generator(long P);
//void primitive_root(long N, long P);
void convolution(long*vec1, long* vec2, long* result, long* temp, long N, long P, long root);
long power_mod(long base, long exp, long mod);
long inverse(long val, long mod);
#endif
