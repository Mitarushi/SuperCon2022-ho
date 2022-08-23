#include <mpi.h>
#include <bits/stdc++.h>
#include "tests/sc_header.h"  // このヘッダーを必ず include すること。

int qs[sc::N_MAX];
int myid, n_procs;
uint8_t table[sc::M_MAX][sc::N_MAX / 8];
short hamming_distance[sc::M_MAX][sc::M_MAX];

constexpr int SIMULATE_BLOCK_SIZE = 128;
void simulate_block(char *s, uint8_t *result, int i_form) {
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
        int index = std::min(i_form + i, sc::N_MAX - 1);
        result[index / 8] |= sc::F[q[i]] << (index % 8);
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
        MPI_Bcast(table + from, sc::N_MAX / 8 * (to - from), MPI_CHAR, id, MPI_COMM_WORLD);
    }
}

void gen_table2() {
    const int bs = (sc::M_MAX + (n_procs - 1)) / n_procs;

    const int index_from = bs * myid;
    const int index_to = std::min(bs * (myid + 1), sc::M_MAX);

#pragma omp parallel for
    for (int i = index_from; i < index_to; i++) {
        for (int j = 0; j < sc::M_MAX; j++) {
            if (i >= j) continue;
            hamming_distance[i][j] = 0;
            for (int k = 0; k < sc::N_MAX / 8; k++) {
                uint8_t t = table[i][k] ^ table[j][k];
                for (int l = 0; l < 8; l++) {
                    hamming_distance[i][j] += (t >> l) & 1;
                }
            }
        }
    }
}

// main関数で入力を読み込んだ後、以下の関数が実行される。
void run() {
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    gen_table();
    gen_table2();

    std::vector<std::tuple<int, int, int>> hamming_sort;
    hammig_sort.reserve(sc::M_MAX * (sc::M_MAX - 1) / 2);

    for (int i = 0; i < sc::M_MAX; i++) {
        for (int j = i + 1; j < sc::M_MAX; j++) {
            hamming_sort.emplace_back(hamming_distance[i][j], i, j);
        }
    }
    std::sort(hamming_sort.begin(), hamming_sort.end());

}

int main(int argc, char **argv) {
    sc::initialize(argc, argv);  // はじめに sc::initialize(argc, argv) を必ず呼び出すこと。

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < sc::N_MAX; j++) {
            sc::T[i][j]--;
        }
    }
    run();

    sc::finalize();  // おわりに sc::finalize() を必ず呼び出すこと。
}