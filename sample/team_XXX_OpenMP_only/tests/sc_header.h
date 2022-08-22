#include <omp.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace sc {
constexpr int CASE_SIZE     = 10;
constexpr double TIME_LIMIT = 60.0;
constexpr int A_SIZE        = 2;

constexpr int N_MIN = 1;
constexpr int N_MAX = 1000;
constexpr int M_MIN = 1;
constexpr int M_MAX = 2000;
constexpr int L_MAX = (M_MAX * 500);

int n;
int T[A_SIZE][N_MAX];
int F[N_MAX];
int m;
char *w[M_MAX];
int w_len[M_MAX];

// インターフェース

// 出力に用いる関数。出力の際は、この関数を必ず使用しなくてはならない。
// int k: 使用する状態の個数を表す。
// int qs[] : 使用する状態の情報を持った長さ k の配列。
void output(int k, int qs[]);

// 現在の実行時間（単位秒）を返す関数。
double get_elapsed_time();

// 実装
#define SC22_KEY 334  // 審査ではSC22_KEYの値を別の「ある値」に変更する。
template <int key> double &time0() {
  static double t;
  return t;
}
double get_elapsed_time() { return omp_get_wtime() - time0<SC22_KEY>(); }

// 入力を読み込む関数。
void input() {
  int rank = 0;
#ifdef MPI_VERSION
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank == 0) {
    if (scanf("%d", &n) != 1) assert(false);
    for (int c = 0; c < A_SIZE; c++) {
      for (int i = 0; i < n; i++) {
        if (scanf("%d", &T[c][i]) != 1) assert(false);
      }
    }
    for (int i = 0; i < n; i++) {
      if (scanf("%d", &F[i]) != 1) assert(false);
    }
    if (scanf("%d", &m) != 1) assert(false);
    char s[L_MAX];
    for (int i = 0; i < m; i++) {
      if (scanf("%s", s) != 1) assert(false);
      w_len[i] = strlen(s);
      w[i]     = (char *)malloc((strlen(s) + 1) * sizeof(char));
      strcpy(w[i], s);
    }
  }
#ifdef MPI_VERSION
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(T, A_SIZE * N_MAX, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(F, N_MAX, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(w_len, M_MAX, MPI_INT, 0, MPI_COMM_WORLD);
  for (int i = 0; i < m; i++) {
    if (rank != 0) MPI_Alloc_mem((w_len[i] + 1) * sizeof(char), MPI_INFO_NULL, &w[i]);
    MPI_Bcast(w[i], w_len[i] + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
#endif
}

// 出力に用いる関数。出力の際は、この関数を必ず使用しなくてはならない。
// int k: 使用する状態の個数を表す。
// int qs[] : 使用する状態の情報を持った長さ k の配列。
void output(int k, int qs[]) {
  int rank = 0;
#ifdef MPI_VERSION
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank) return;
  if (omp_get_thread_num()) return;

  printf("%d : begin_output : %f\n", SC22_KEY, get_elapsed_time());
  printf("%d : %d\n", SC22_KEY, k);
  printf("%d : ", SC22_KEY);
  for (int j = 0; j < k; j++) printf("%d ", qs[j]);
  printf("\n");
  printf("%d : end_output : %f\n", SC22_KEY, get_elapsed_time());
  fflush(stdout);
}

void initialize(int &argc, char **&argv) {
#ifdef MPI_VERSION
  MPI_Init(&argc, &argv);
#endif
  input();
  time0<SC22_KEY>() = omp_get_wtime();
}

void finalize() {
  for (int i = 0; i < m; i++) { free(w[i]); }
#ifdef MPI_VERSION
  MPI_Finalize();
#endif
}

#undef SC22_KEY

}  // namespace sc
