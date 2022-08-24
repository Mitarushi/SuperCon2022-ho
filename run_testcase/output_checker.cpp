// set PATH=C:\msys64\mingw64\bin;C:\msys64\usr\local\bin;C:\msys64\usr\bin;C:\msys64\bin;%PATH%
// set MSYSTEM=MINGW64
// bash
//g++ 6.cpp -std=c++2a -O2 -I .
#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using ld = long double;
 
using vi = vector<int>;
using vvi = vector<vi>;
using vll = vector<ll>;
using vvll = vector<vll>;
using vld = vector<ld>;
using vvld = vector<vld>;
using vst = vector<string>;
using vvst = vector<vst>;
 
#define fi first
#define se second
#define pb push_back
#define eb emplace_back
#define pq_big(T) priority_queue<T,vector<T>,less<T>>
#define pq_small(T) priority_queue<T,vector<T>,greater<T>>
#define all(a) a.begin(),a.end()
#define rep(i,start,end) for(ll i=start;i<(ll)(end);i++)
#define per(i,start,end) for(ll i=start;i>=(ll)(end);i--)
#define uniq(a) sort(all(a));a.erase(unique(all(a)),a.end())

constexpr unsigned long long SEED = 1ull;

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
unsigned long long range_rnd(unsigned long long l, unsigned long long r)
{ return l + rnd(r - l + 1); }


bool output_check(int n, int m, std::vector<int> &ans, std::vector<std::vector<int>> &accept, std::vector<std::vector<int>> &hamming_distance){
    for (int i = 0; i < m; i++) {
        for (int j = i + 1; j < m; j++) {
            if (hamming_distance[i][j] == 0) {
                continue;
            }

            bool equiv = 1;
            for (int index : ans) {
                if (accept[i][index] != accept[j][index]) {
                    equiv = 0;
                    break;
                }
            }

            if (equiv == 0) {
                continue;
            }

            for (int k = 0; k < n; k++) {
                if (accept[i][k] != accept[j][k]) {
                    // cout << "ERROR : "<< i << " " << j << endl;
                    return 0;
                }
            }
        }
    }
    return 1;
}