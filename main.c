#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

size_t reverse_bits(int log2N, size_t input) {
  if (log2N > sizeof(size_t) * CHAR_BIT) {
    fprintf(stderr, "log2N = %d too large\n", log2N);
    exit(1);
  }
  size_t ans = 0;
  for (int i = 0; i < log2N; i++) {
    if (input & (1 << i)) {
      ans |= 1 << (log2N - i - 1);
    }
  }
  return ans;
}

void butterfly(double complex *pa, double complex *pb, double complex w) {
    double complex new_a = *pa + *pb;
    double complex new_b = (*pa - *pb) * w;
    *pa = new_a;
    *pb = new_b;
}

double complex W(int m, int N) { return cexp(-2 * atan(1) * 4 * I * m / N); }

void fft(int N, const double complex in_[N], double complex out[N]) {
  double complex in[N];
  for (int i = 0; i < N; i++) {
    in[i] = in_[i];
  }
  int log2_N = log2(N);
  if (exp2(log2_N) != N) {
    fprintf(stderr, "Not a power of 2: %d\n", N);
    exit(1);
  }

  for (int stage = 0; stage < log2_N; stage++) {
    int mask = N >> (stage + 1);
    for (int i = 0; i < N; i++) {
      if (i & mask) continue;
      butterfly(&in[i], &in[i | mask], W(i & (mask - 1), N >> stage));
    }
  }

  for (int i = 0; i < N; i++) {
    out[i] = in[reverse_bits(log2(N), i)];
  }
}

// out_k = \sum_{m=0}^{N-1} in_m \exp(-2\pi imk / N)
void dft(int N, const double complex in[N], double complex out[N]) {
  const double PI = atan(1) * 4;
  for (int k = 0; k < N; k++) {
    out[k] = 0.0;
    for (int m = 0; m < N; m++) {
      out[k] += in[m] * cexp(-2 * PI * I * m * k / N);
    }
  }
}

int main(void) {
  int N = 16;
  const double complex in[16] = {3, -1, 4, 1, -5, 9, 2, 6, -5, 3, 5+8 *I, 9, -7, 9, 3, 2};
  double complex out[N];
  dft(N, in, out);
  for (int i = 0; i < N; i++) {
    printf("dft: out[%d] = %f + %fi\n", i, creal(out[i]), cimag(out[i]));
  }
  printf("\n");

  fft(N, in, out);
  for (int i = 0; i < N; i++) {
    printf("fft: out[%d] = %f + %fi\n", i, creal(out[i]), cimag(out[i]));
  }
}