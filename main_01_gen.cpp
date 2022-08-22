#include <mpi.h>

#include "tests/sc_header.h"  // このヘッダーを必ず include すること。

int qs[sc::N_MAX];
int myid, n_procs;
bool table[sc::M_MAX][sc::N_MAX];

constexpr int SIMULATE_BLOCK_SIZE = 32;
void simulate_block(char *s, bool *result, int i_form) {
    short q[SIMULATE_BLOCK_SIZE];
    for (int i = 0; i < SIMULATE_BLOCK_SIZE; i++) {
        q[i] = std::min(i_form + i, sc::N_MAX - 1);
    }

    for (int k = 0; s[k] != '\0'; k++) {
        int c = s[k] - 'a';

        for (int i = 0; i < SIMULATE_BLOCK_SIZE; i++) {
            q[i] = sc::T[c][q[i]];
        }
    }

    for (int i = 0; i < SIMULATE_BLOCK_SIZE; i++) {
        result[std::min(i_form + i, sc::N_MAX - 1)] = sc::F[q[i]];
    }
}

void gen_table() {
    const int bs = (sc::M_MAX + (n_procs - 1)) / n_procs;

    const int index_from = bs * myid;
    const int index_to = std::min(bs * (myid + 1), sc::M_MAX);

#pragma omp parallel for
    for (int i = index_from; i < index_to; i++) {
        for (int j = 0; j < sc::N_MAX; j += SIMULATE_BLOCK_SIZE) {
            simulate_block(sc::w[i], table[i], j);
        }
    }

    for (int id = 0; id < n_procs; id++) {
        int from = bs * id;
        int to = std::min(bs * (id + 1), sc::M_MAX);
        MPI_Bcast(table + from, sc::N_MAX * (to - from), MPI_C_BOOL, id,
                  MPI_COMM_WORLD);

        int cnt = 0;
        for (int i = from; i < to; i++) {
            for (int j = 0; j < sc::N_MAX; j++) {
                cnt += table[i][j];
            }
        }
    }
}

// main関数で入力を読み込んだ後、以下の関数が実行される。
void run() {
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    gen_table();

    int cnt = 0;
    for (int i = 0; i < sc::M_MAX; i++) {
        for (int j = 0; j < sc::N_MAX; j++) {
            cnt += table[i][j];
        }
    }

    printf("%d\n", cnt);
}

int main(int argc, char **argv) {
    sc::initialize(
        argc,
        argv);  // はじめに sc::initialize(argc, argv) を必ず呼び出すこと。

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < sc::N_MAX; j++) {
            sc::T[i][j]--;
        }
    }
    run();

    sc::finalize();  // おわりに sc::finalize() を必ず呼び出すこと。
}