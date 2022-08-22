#include <cstdio>
#include <cstdlib>
// 審査ではSEEDの値を別の「ある値」に変更する。
constexpr unsigned long long SEED = 1ull;

constexpr int CASE_SIZE = 10;
constexpr int A_SIZE    = 2;
constexpr int N         = 1000;
constexpr int M         = 2000;
constexpr int L         = (M * 500);

unsigned long long xor_shift() {
  static unsigned long long x = SEED;
  x                           = x ^ (x << 13);
  x                           = x ^ (x >> 7);
  x                           = x ^ (x << 17);
  return x;
}

// [0, x)
unsigned long long rnd(unsigned long long x) { return xor_shift() % x; }

// [l, r]
unsigned long long range_rnd(unsigned long long l, unsigned long long r) { return l + rnd(r - l + 1); }

int main() {
  for (int t = 1; t <= CASE_SIZE; t++) {
    char file_name[64];
    sprintf(file_name, "random_%02d.in", t);

    int **T = (int **)malloc(A_SIZE * sizeof(int *));
    for (int c = 0; c < A_SIZE; c++) {
      T[c] = (int *)malloc(N * sizeof(int));
      for (int i = 0; i < N; i++) { T[c][i] = range_rnd(1, N); }
    }
    int *F = (int *)malloc(N * sizeof(int));
    for (int i = 0; i < N; i++) { F[i] = range_rnd(0, 1); }

    int *ls = (int *)calloc(M, sizeof(int));
    for (int i = 0; i < L; i++) {
      int pos = i < M ? i : range_rnd(0, M - 1);
      ls[pos]++;
    }
    char **ws = (char **)malloc(M * sizeof(char *));
    for (int i = 0; i < M; i++) {
      ws[i] = (char *)malloc((ls[i] + 1) * sizeof(char));
      for (int j = 0; j < ls[i]; j++) { ws[i][j] = 'a' + range_rnd(0, A_SIZE - 1); }
      ws[i][ls[i]] = '\0';
    }

    {
      FILE *fp = fopen(file_name, "w");
      if (fp == NULL) return 1;
      {
        fprintf(fp, "%d\n", N);
        for (int c = 0; c < A_SIZE; c++) {
          for (int i = 0; i < N; i++) {
            if (i > 0) fprintf(fp, " ");
            fprintf(fp, "%d", T[c][i]);
          }
          fprintf(fp, "\n");
        }
        for (int i = 0; i < N; i++) {
          if (i > 0) fprintf(fp, " ");
          fprintf(fp, "%d", F[i]);
        }
        fprintf(fp, "\n");
      }
      {
        fprintf(fp, "%d\n", M);
        for (int i = 0; i < M; i++) { fprintf(fp, "%s\n", ws[i]); }
      }
      fclose(fp);
    }

    free(ls);
    for (int i = 0; i < M; i++) free(ws[i]);
    free(ws);
    free(F);
    for (int c = 0; c < A_SIZE; c++) free(T[c]);
    free(T);
  }
  return 0;
}