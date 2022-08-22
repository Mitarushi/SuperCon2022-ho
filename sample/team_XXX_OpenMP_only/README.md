# README

OpenMPのみを使用する場合のサンプル。

提出ファイルは `run_cases.sh` と `main.cpp` の２つのファイル。

`run_cases.sh` のパラメータおよび `main.cpp` の中身のみを書き換え、これら以外（`tests` 以下）を書き換えないで正しく動作するものを提出すること。

## 構成

    main.cpp         # OpenMPのみを使用する場合のサンプルプログラム（提出ファイル）

    run_cases.sh     # ジョブを投げるためのシェルスクリプト（提出ファイル）

    README.md        # 説明のためのファイル（このファイル）

    tests            # 提出ファイル以外のファイルが入ったディレクトリ（このディレクトリ以下のファイルを書き換えても審査に反映されないため注意すること）

    log              # tests/jobscript.sh を実行すると出力結果のファイルが入るディレクトリ

tests 以下

    sc_header.h      # ヘッダファイル

    jobscript.sh     # main.cpp をコンパイルして、tests/random_01.in から tests/random_10.in の入力例を実行するスクリプト

    00_sample_01.in  # 問題文の入力例1

    00_sample_02.in  # 問題文の入力例2

    generator.cpp    # 入力生成器のプログラム

    judge.cpp        # ジャッジのためのプログラム

## ジャッジスクリプトの実行方法（富岳の使用方法についてはチュートリアルの各資料を参考にすること）

1. シェルスクリプトに実行権限を付与

        chmod +x run_cases.sh tests/jobscript.sh

2. pjsub でジョブスクリプト `run_cases.sh` を投げる

    以下のように実行する（`run_cases.sh` 内の各パラメータを書き換えて使用する）。

        pjsub -L rscgrp=excl_supercon_2208-1,freq=2200,eco_state=2 run_cases.sh

3. 実行後、`log/result.${PJM_JOBID}.out` に各ケースのスコアと実行時間が出力される（${PJM_JOBID}は、ジョブ毎に与えられるある整数値）。
