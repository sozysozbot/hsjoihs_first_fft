#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


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
  const double complex in[8] =
  { 3,
    -1,
    4,
    1,
    -5,
    9,
    2,
    6 };

  double complex out[8];
  dft(8, in, out);

  for (int i = 0; i < 8; i++) {
    printf("out[%d] = %f + %fi\n", i, creal(out[i]), cimag(out[i]));
  }
}