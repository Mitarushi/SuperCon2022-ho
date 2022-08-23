#include <bits/stdc++.h>
#include <mpi.h>

#include "tests/sc_header.h"  // このヘッダーを必ず include すること。

int qs[sc::N_MAX];
int myid, n_procs;
bool table[sc::M_MAX][sc::N_MAX];
uint16_t hamming_distance[sc::M_MAX][sc::M_MAX];

constexpr int SIMULATE_BLOCK_SIZE = 128;
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
        int index = std::min(i_form + i, sc::N_MAX - 1);
        result[index] = sc::F[q[i]];
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
        MPI_Bcast(table + from, sc::N_MAX * (to - from), MPI_C_BOOL, id, MPI_COMM_WORLD);
    }
}

void gen_table2() {
    const int bs = (sc::M_MAX + (n_procs - 1)) / n_procs;

    const int index_from = bs * myid;
    const int index_to = std::min(bs * (myid + 1), sc::M_MAX);

#pragma omp parallel for
    for (int i = index_from; i < index_to; i++) {
        for (int j = 0; j < sc::M_MAX; j++) {
            hamming_distance[i][j] = 0;
            for (int k = 0; k < sc::N_MAX; k++) {
                hamming_distance[i][j] += table[i][k] ^ table[j][k];
            }
        }
    }

    for (int id = 0; id < n_procs; id++) {
        int from = bs * id;
        int to = std::min(bs * (id + 1), sc::M_MAX);
        MPI_Bcast(hamming_distance + from, sc::M_MAX * (to - from), MPI_SHORT, id, MPI_COMM_WORLD);
    }
}

struct hash {
    uint64_t x;

    hash() : x(0) {}

    void next(int c) {
        x = x ^ (x << 33);
        x = x ^ (x >> 27);
        x = x ^ (x << 5);
        x ^= c;
    }
};

std::vector<std::vector<int>> get_small_hamming(int k) {
    std::vector<std::tuple<int, int, int>> hamming_sort;
    hamming_sort.reserve(sc::M_MAX * (sc::M_MAX - 1) / 2);

    for (int i = 0; i < sc::M_MAX; i++) {
        for (int j = i + 1; j < sc::M_MAX; j++) {
            if (hamming_distance[i][j] == 0) {
                continue;
            }
            hamming_sort.emplace_back(hamming_distance[i][j], i, j);
        }
    }
    std::sort(hamming_sort.begin(), hamming_sort.end());

    std::vector<std::vector<int>> result;
    result.reserve(k);

    std::unordered_set<uint64_t> hash_set;
    for (auto [c, p, q] : hamming_sort) {
        hash h;
        std::vector<int> r;
        for (int i = 0; i < sc::N_MAX; i++) {
            int t = table[p][i] ^ table[q][i];
            if (t) {
                h.next(i);
                r.push_back(i);
            }
        }

        if (!hash_set.count(h.x)) {
            hash_set.insert(h.x);
            result.push_back(r);
        }

        if (result.size() == k) break;
    }

    return result;
}

// main関数で入力を読み込んだ後、以下の関数が実行される。
void run() {
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    gen_table();
    gen_table2();

    const int k = 1000;
    std::vector<std::vector<int>> small_hamming = get_small_hamming(k);

    std::vector<int> greedy;
    std::vector<char> satisfied(k, 0);
    while (true) {
        std::vector<double> weight_count(sc::N_MAX, 0.0);

        for (int i = 0; i < k; i++) {
            if (satisfied[i]) {
                continue;
            }

            double inv = 1.0 / small_hamming[i].size();
            for (int j : small_hamming[i]) {
                weight_count[j] += inv;
            }
        }

        int it = std::max_element(weight_count.begin(), weight_count.end()) - weight_count.begin();
        if (weight_count[it] == 0.0) {
            break;
        }

        for (int i = 0; i < k; i++) {
            if (satisfied[i]) {
                continue;
            }

            if (std::count(small_hamming[i].begin(), small_hamming[i].end(), it)) {
                satisfied[i] = 1;
            }
        }

        greedy.push_back(it);
    }

    while (true) {
        int cnt = 0;

        std::vector<double> weight_count(sc::N_MAX, 0.0);

        for (int i = 0; i < sc::M_MAX; i++) {
            for (int j = 0; j < i; j++) {
                if (hamming_distance[i][j] == 0) {
                    continue;
                }

                // std::cout << i << " " << j << std::endl;

                bool flag = true;
                for (int k : greedy) {
                    if (table[i][k] ^ table[j][k]) {
                        flag = false;
                        break;
                    }
                }
                if (!flag) {
                    continue;
                }

                cnt += 1;

                std::vector<int> use;
                for (int k = 0; k < sc::N_MAX; k++) {
                    if (table[i][k] ^ table[j][k]) {
                        use.push_back(k);
                    }
                }

                double inv = 1.0 / use.size();
                for (int k : use) {
                    weight_count[k] += inv;
                }
            }
        }

        if (cnt == 0) {
            break;
        }

        int it = std::max_element(weight_count.begin(), weight_count.end()) - weight_count.begin();
        greedy.push_back(it);
    }

    if (myid == 0) {
        sc::output(greedy.size(), greedy.data());
    }
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