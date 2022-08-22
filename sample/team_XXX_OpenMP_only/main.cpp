// #include <mpi.h>  // MPI を使用する場合、このヘッダーを必ず include すること。

#include "tests/sc_header.h"  // このヘッダーを必ず include すること。

int qs[sc::N_MAX];

// main関数で入力を読み込んだ後、以下の関数が実行される。
void run() {
  for (int i = 0; i < sc::n; i++) qs[i] = i + 1;
  sc::output(sc::n, qs);
}

int main(int argc, char** argv) {
  sc::initialize(argc, argv);  // はじめに sc::initialize(argc, argv) を必ず呼び出すこと。
  run();
  sc::finalize();  // おわりに sc::finalize() を必ず呼び出すこと。
}
