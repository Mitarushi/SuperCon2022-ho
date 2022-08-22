// ジャッジプログラムのサンプル。
#include <bits/stdc++.h>
using namespace std;

#include "sc_header.h"
using namespace sc;

#define SC22_KEY 334

// BEGIN: This code block is copied from ../main.cpp

// 疑似乱数の生成に使用する（これ以外の疑似乱数を用いてもよい）。
struct xor_shift {
  unsigned long long x;
  xor_shift(unsigned long long SEED) { x = SEED; }
  unsigned long long next() {
    x = x ^ (x << 13);
    x = x ^ (x >> 7);
    x = x ^ (x << 17);
    return x;
  }
  // [0, x)
  unsigned long long rnd(unsigned long long x) { return next() % x; }
  // [l, r]
  unsigned long long range_rnd(unsigned long long l, unsigned long long r) { return l + rnd(r - l + 1); }
};

// accs[j][i] = 状態 i からはじめて文字列 w_j の遷移関数で遷移した先の状態が受理状態の場合 1、受理状態でない場合 0
// が代入される。
int accs[sc::M_MAX][sc::N_MAX];
// オートマトンの遷移関数をシミュレーションする関数
void simulate() {
  for (int j = 0; j < sc::m; j++) {
    int i, q, k;
#pragma omp parallel for private(q, k)
    for (i = 0; i < sc::n; i++) {
      q = i;
      for (k = 0; sc::w[j][k] != '\0'; k++) { q = sc::T[sc::w[j][k] - 'a'][q]; }
      accs[j][i] = sc::F[q];
    }
  }
}

// 状態 qs[0] から 状態 qs[k-1] のk個の状態を使った時に区別できる文字列のペアの個数を返す関数。
int get_sep_c(int k, int *qs) {
  int c = 0;
  for (int j1 = 0; j1 < sc::m; j1++) {
    for (int j2 = j1; j2 < sc::m; j2++) {
      bool ok = false;
      for (int i = 0; i < k; i++) {
        if (accs[j1][qs[i] - 1] != accs[j2][qs[i] - 1]) {
          ok = true;
          break;
        }
      }
      c += ok;
    }
  }
  return c;
}
// 全ての状態を使用する場合に区別できる文字列のペアの個数を返す関数。
int get_max_sep_c() {
  int qs[sc::N_MAX];
  for (int i = 0; i < sc::n; i++) qs[i] = 1 + i;
  return get_sep_c(sc::n, qs);
}

// スコアを計算する関数（ただし、出力が正しいかどうかのチェックはここではしない）。
// int k: 使用する状態の個数を表す。
// int *qs: 使用する状態の情報を持った長さ k の配列。
long get_score(int k, int *qs) {
  long score = k * sc::n * sc::n;
  for (int i = 0; i < k; i++) score += qs[i];
  return score;
}
long get_worst_score() { return 2ll * sc::n * sc::n * sc::n; }

// END: This code block is copied from ../main.cpp

double ouf_time = sc::TIME_LIMIT;
int ouf_k       = -1;
int ouf_qs[N_MAX];
void read_out(const char *file) {
  ifstream ifs(file);
  if (!ifs) return;
  while (true) {
    int tmp_k;
    int tmp_qs[N_MAX];
    string _delim1, _delim2, _delim3, _delim4;
    double tmp_t_begin, tmp_t_end;
    // line 1
    {
      ifs >> _delim1 >> _delim2 >> _delim3 >> _delim4;
      if (_delim1 == "") break;
      if (_delim1 + _delim2 != to_string(SC22_KEY) + ":") break;
      if (_delim3 + _delim4 != "begin_output:") break;
      ifs >> tmp_t_begin;
    }
    // line 2
    {
      ifs >> _delim1 >> _delim2 >> tmp_k;
      if (_delim1 + _delim2 != to_string(SC22_KEY) + ":") break;
      if (!(1 <= tmp_k && tmp_k <= n)) break;
    }
    // line 3
    {
      ifs >> _delim1 >> _delim2;
      if (_delim1 + _delim2 != to_string(SC22_KEY) + ":") break;
      for (int i = 0; i < tmp_k; i++) {
        ifs >> tmp_qs[i];
        if (!(1 <= tmp_qs[i] && tmp_qs[i] <= n)) break;
      }
    }
    // line 4
    {
      ifs >> _delim1 >> _delim2 >> _delim3 >> _delim4;
      if (_delim1 + _delim2 != to_string(SC22_KEY) + ":") break;
      if (_delim3 + _delim4 != "end_output:") break;
      if (ifs.eof()) break;
      ifs >> tmp_t_end;
      if (tmp_t_end > sc::TIME_LIMIT) break;
    }

    // update
    ouf_k = tmp_k;
    for (int i = 0; i < tmp_k; i++) ouf_qs[i] = tmp_qs[i];
    ouf_time = tmp_t_end;
  }
}

void judge() {
  const int max_sep_c = get_max_sep_c();
  const int sep_c     = get_sep_c(ouf_k, ouf_qs);
  if (sep_c == max_sep_c) {
    cout << get_score(ouf_k, ouf_qs) << endl;
    cout << ouf_time << endl;
  } else {
    cout << get_worst_score() << endl;
    cout << ouf_time << endl;
  }
}

void run() {
  for (int c = 0; c < sc::A_SIZE; c++)
    for (int i = 0; i < sc::n; i++) sc::T[c][i]--;
  simulate();
  judge();
}

bool is_valid() {
  if (!(1 <= ouf_k && ouf_k <= n)) return false;
  for (int i = 0; i < ouf_k; i++) {
    if (!(1 <= ouf_qs[i] && ouf_qs[i] <= n)) return false;
  }
  return true;
}

int main(int argc, char **argv) {
  assert(argc == 2);
  sc::initialize(argc, argv);
  read_out(argv[1]);
  if (!is_valid()) {
    cout << get_worst_score() << endl;
    cout << ouf_time << endl;
  } else {
    run();
  }
  sc::finalize();
}
