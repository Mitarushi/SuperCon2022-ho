#include <bits/stdc++.h>
#include <mpi.h>

#include "tests/sc_header.h"  // このヘッダーを必ず include すること。

int qs[sc::N_MAX];
int myid, n_procs;
bool table[sc::M_MAX][sc::N_MAX];
std::uint64_t t_table_bit[sc::N_MAX][sc::M_MAX / 64 + 1] = {};
uint16_t hamming_distance[sc::M_MAX][sc::M_MAX];
std::vector<std::pair<int, int>> small_hamming;
int leader[sc::M_MAX];
int leader_size = 0;

const int N_SQUARE = sc::N_MAX * sc::N_MAX;

unsigned long long SEED = 1ull;
unsigned long long xor_shift() {
    static unsigned long long x = SEED;
    x = x ^ (x << 13);
    x = x ^ (x >> 7);
    x = x ^ (x << 17);
    return x;
}

// [0, x)
unsigned long long rnd(unsigned long long x) { return xor_shift() % x; }
// [l, r]
unsigned long long range_rnd(unsigned long long l, unsigned long long r) { return l + rnd(r - l + 1); }

// [0, x)
template <unsigned long long x>
unsigned long long rnd() {
    return xor_shift() % x;
}

void shuffle(std::vector<int> &v) {
    for (int i = 0; i < (int)v.size(); i++) {
        int j = rnd(i + 1);
        std::swap(v[i], v[j]);
    }
}

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

    leader_size = 0;
    for (int i = 0; i < sc::M_MAX; i++) {
        bool t = true;
        for (int j = 0; j < i; j++) {
            if (hamming_distance[i][j] == 0) {
                t = false;
                break;
            }
        }

        if (t) {
            leader[leader_size++] = i;
        }
    }
}

void gen_table3() {
    for (int i = 0; i < leader_size; i++) {
        for (int j = 0; j < sc::N_MAX; j++) {
            t_table_bit[j][i / 64] |= (std::uint64_t(table[leader[i]][j]) << (i % 64));
        }
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
        x *= 35367321;
    }
};

void get_small_hamming(int k) {
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

    // small_hamming;

    std::unordered_set<uint64_t> hash_set;
    for (auto [c, p, q] : hamming_sort) {
        if (c > k) {
            break;
        }
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
            small_hamming.emplace_back(p, q);
        }
    }

    // std::cout << "small_hamming size: " << small_hamming.size() << std::endl;
}

// ハミング距離が短い上位 k ペア採用
void get_small_hamming_for_ranking(int k) {
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

    // small_hamming;

    std::unordered_set<uint64_t> hash_set;
    for (auto [c, p, q] : hamming_sort) {
        if (small_hamming.size() >= k) {
            break;
        }
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
            small_hamming.emplace_back(p, q);
        }
    }

    // std::cout << "small_hamming size: " << small_hamming.size() << std::endl;
}

struct timer {
    double start;
    double end;
    double acc;
    timer() : start(0), end(0), acc(0) {}
    void start_timer() { start = sc::get_elapsed_time(); }
    void end_timer() {
        end = sc::get_elapsed_time();
        acc += end - start;
    }
    double get_time() { return acc; }
};

inline bool is_distinguishable(int i, int j, std::vector<int> &x) {
    if (x.size() > 0) {
        int k = x[x.size() - 1];
        if (table[i][k] ^ table[j][k]) {
            return true;
        }
    }
    for (int k : x) {
        if (table[i][k] ^ table[j][k]) {
            return true;
        }
    }
    return false;
}

timer ti;
std::vector<int> construct(std::vector<int> &from, int priority, int restriction) {
    // std::cout << 0 << std::endl;
    std::vector<int> result = from;

    const int buf_size = 11000;
    std::pair<int, int> buffer[buf_size] = {};
    int buf_idx = 0;
    bool buf_useable = false;

    while (true) {
        std::vector<short> count(sc::N_MAX, 0);

        bool is_end = true;
        bool is_small = result.size() < 25;

        if (is_small) {
            // std::cout << small_hamming.size() << std::endl;
            int sh_size = small_hamming.size();
#pragma omp parallel for
            for (int p = 0; p < sh_size; p++) {
                int i = small_hamming[p].first;
                int j = small_hamming[p].second;
                bool distinguish_able = is_distinguishable(i, j, result);

                if (distinguish_able) {
                    continue;
                }

                is_end = false;

                for (int k = 0; k < sc::N_MAX; k++) {
                    if (table[i][k] ^ table[j][k]) {
#pragma omp atomic
                        count[k]++;
                    }
                }
            }
        }
        if (is_end || !is_small) {
            if (buf_useable) {
                int idx = 0;
                int end = buf_idx - 1;
                for (int p = 0; p < buf_idx; p++) {
                    int i = buffer[idx].first;
                    int j = buffer[idx].second;
                    if (!is_distinguishable(i, j, result)) {
                        idx++;
                        continue;
                    }

                    std::swap(buffer[idx], buffer[end]);
                    end--;
                }

                buf_idx = idx;

            } else {
                buf_idx = 0;
                ti.start_timer();
#pragma omp parallel for schedule(dynamic)
                for (int i_idx = 0; i_idx < leader_size; i_idx++) {
                    int i = leader[i_idx];

                    std::uint64_t is_distinguishable[sc::M_MAX / 64 + 1] = {};

                    for (int k : result) {
                        std::uint64_t mask = 0ULL - table[i][k];
                        for (int j_idx = 0; j_idx < i_idx / 64 + 1; j_idx++) {
                            is_distinguishable[j_idx] |= mask ^ t_table_bit[k][j_idx];
                        }
                    }

                    for (int j_idx = 0; j_idx < i_idx; j_idx++) {
                        if (!((is_distinguishable[j_idx / 64] >> (j_idx % 64)) & 1)) {
                            int j = leader[j_idx];

                            int t;

#pragma omp atomic capture
                            t = buf_idx++;
                            // t = buf_idx - 1;

                            if (t < buf_size) {
                                buffer[t] = std::make_pair(i, j);
                            }
                        }
                    }
                }

                buf_useable = buf_idx < buf_size;
                ti.end_timer();
            }

#pragma omp parallel for
            for (int idx = 0; idx < buf_idx; idx++) {
                int i = buffer[idx].first;
                int j = buffer[idx].second;

                is_end = false;

                for (int k = 0; k < sc::N_MAX; k++) {
                    if (table[i][k] ^ table[j][k]) {
#pragma omp atomic
                        count[k]++;
                    }
                }
            }
        }

        if (is_end) {
            break;
        }

        std::vector<double> ratio(sc::N_MAX);
        for (int i = 0; i < sc::N_MAX; i++) {
            ratio[i] = double(N_SQUARE + i + 1) / std::pow(count[i] + 1e-9, 1.5);
        }

        int min_index = std::min_element(ratio.begin(), ratio.end()) - ratio.begin();

        int p = rnd<100>();
        if (p > priority) {
            std::vector<int> ok;
            double upper = ratio[min_index] * (1 + (double)restriction / 100.0);
            for (int i = 0; i < sc::N_MAX; i++) {
                if (ratio[i] < upper) {
                    ok.push_back(i);
                }
            }

            // std::cout << "ok size: " << ok.size() << std::endl;

            min_index = ok[rnd(ok.size())];
        }

        // std::cout << result.size() << " " << min_index << " " << ratio[min_index] << std::endl;
        result.push_back(min_index);
    }

    // std::cout << "construct time: " << t.get_time() << std::endl;

    // std::cout << result.size() << std::endl;
    // std::cout << 4 << std::endl;

    return result;
}

std::int64_t get_score(std::vector<int> &x) {
    std::int64_t s = 0;
    for (int i : x) {
        s += N_SQUARE + i + 1;
    }
    return s;
}

std::vector<int> neighbor_search(std::vector<int> &from, int priority, int restriction, int magnitude) {
    // std::cout << 1 << std::endl;
    std::vector<int> result = from;

    for (int iter = 0; iter < 100; iter++) {
        std::int64_t prev_score = get_score(result);
        std::vector<int> prev_result = result;

        int remove_num = (result.size() * magnitude + 99) / 100;
        shuffle(result);
        result.resize(result.size() - remove_num);

        std::vector<int> result2 = construct(result, priority, restriction);
        result.insert(result.end(), result2.begin(), result2.end());

        // std::cout << result.size() << " " << result2.size() << std::endl;

        std::int64_t score = get_score(result);

        if (prev_score < score) {
            result = prev_result;
        }
    }

    // std::cout << 6 << std::endl;

    return result;
}

inline int linear_inter(int start, int end) {
    double t = sc::get_elapsed_time() / sc::TIME_LIMIT;
    t = std::min(1.0, t);
    return start + std::round((end - start) * t);
}

std::int64_t optimal_score = 1000LL * N_SQUARE;

std::vector<int> set_cover(std::vector<int> &prev_best) {
    // std::cout << 2 << std::endl;
    const int priority = linear_inter(15, 30);
    const int restriction = linear_inter(10, 10);
    const int magnitude = linear_inter(20, 19);
    const int prev_remain = linear_inter(65, 80);
    const int prev_ratio = linear_inter(40, 80);

    double upper = 1.0 + (double)restriction / 100.0;

    std::vector<int> result;

    // std::cout << "start" << std::endl;

    for (int iter = 0; iter < 3; iter++) {
        std::vector<int> empty;

        if (rnd<100>() < prev_ratio) {
            empty = prev_best;
            shuffle(empty);
            empty.resize(empty.size() * prev_remain / 100);
        }

        std::vector<int> x = construct(empty, priority, restriction);
        std::int64_t score = get_score(x);

        // std::cout << score << " " << optimal_score << std::endl;

        if (score < optimal_score) {
            optimal_score = score;
            result = x;
        }

        if (score < optimal_score * upper) {
            x = neighbor_search(x, priority, restriction, magnitude);

            score = get_score(x);

            // std::cout << "!!   " << score << " " << optimal_score << std::endl;
            if (score < optimal_score) {
                optimal_score = score;
                result = x;
            }
        }
    }
    // std::cout << 5 << std::endl;

    return result;
}

bool output_check(int n, int m, std::vector<int> &ans) {
    for (int i = 0; i < m; i++) {
        for (int j = i + 1; j < m; j++) {
            if (hamming_distance[i][j] == 0) {
                continue;
            }

            bool equiv = 1;
            for (int index : ans) {
                if (table[i][index] != table[j][index]) {
                    equiv = 0;
                    break;
                }
            }

            if (equiv == 0) {
                continue;
            }

            for (int k = 0; k < n; k++) {
                if (table[i][k] != table[j][k]) {
                    // cout << "ERROR : "<< i << " " << j << endl;
                    return 0;
                }
            }
        }
    }
    return 1;
}

// main関数で入力を読み込んだ後、以下の関数が実行される。
void run() {
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // std::cout << "pid" << getpid() << " " << myid << std::endl;

    SEED = myid * 314159265 + 358979323;

    gen_table();
    gen_table2();
    gen_table3();
    get_small_hamming(55);

    std::vector<int> best;
    int best_array[sc::N_MAX];

    int cnt = 0;

    while (sc::get_elapsed_time() < sc::TIME_LIMIT) {
        long prev_optimal_score = optimal_score;
        std::vector<int> result = set_cover(best);
        // std::cout << 'p' << myid << " " << optimal_score << std::endl;
        // std::cout << "acc: " << ti.acc << " " << optimal_score << " cnt " << cnt++ << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);

        struct score_id {
            long score;
            int id;
        } cur, mv;
        cur = {optimal_score, myid};

        if (result.size() == 0) {
            cur.score++;
        }

        MPI_Allreduce(&cur, &mv, 1, MPI_LONG_INT, MPI_MINLOC, MPI_COMM_WORLD);

        if (prev_optimal_score < mv.score) {
            continue;
        }
        optimal_score = mv.score;

        if (mv.id == myid) {
            best = result;
            for (int i = 0; i < best.size(); i++) {
                best_array[i] = best[i];
            }
        }

        MPI_Bcast(best_array, sc::N_MAX, MPI_INT, mv.id, MPI_COMM_WORLD);

        best.resize(mv.score / N_SQUARE);
        for (int i = 0; i < best.size(); i++) {
            best[i] = best_array[i];
        }

        // std::cout << 'q' << myid << " " << optimal_score << std::endl;

        if (myid == 0) {
            for (int &i : best) {
                i++;
            }
            sc::output(best.size(), best.data());
            for (int &i : best) {
                i--;
            }

            bool is_ok = output_check(sc::N_MAX, sc::M_MAX, best);
            if (is_ok) {
                // std::cout << "OK" << std::endl;
            } else {
                // std::cout << "ERROR" << std::endl;
            }
        }
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