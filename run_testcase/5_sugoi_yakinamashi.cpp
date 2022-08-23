// g++ 5.cpp -std=c++17 -O2 -I .
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
#define pq_big(T) priority_queue<T, vector<T>, less<T>>
#define pq_small(T) priority_queue<T, vector<T>, greater<T>>
#define all(a) a.begin(), a.end()
#define rep(i, start, end) for (ll i = start; i < (ll)(end); i++)
#define per(i, start, end) for (ll i = start; i >= (ll)(end); i--)
#define uniq(a)   \
    sort(all(a)); \
    a.erase(unique(all(a)), a.end())

constexpr unsigned long long SEED = 1ull; // 1ull,2ull なら 39 / 5ull なら 38

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

void shuffle(vector<int> &v) {
    for (int i = 0; i < (int)v.size(); i++) {
        int j = rnd(i + 1);
        swap(v[i], v[j]);
    }
}

std::vector<int> construction(std::vector<int> &uncover, std::vector<std::vector<char>> &count, int priority, int restriction) {
    std::vector<int> result;
    std::vector<int> uncover2 = uncover;

    while (uncover2.size() > 0) {
        std::vector<double> ratio(1000, 1e9);
        for (int i = 0; i < 1000; i++) {
            double k = 1e-9;
            for (int j : uncover2) {
                k += count[i][j];
            }
            ratio[i] = 1 / k;
        }

        int min_index = std::min_element(ratio.begin(), ratio.end()) - ratio.begin();

        int p = rnd(100);
        if (p > priority) {
            std::vector<int> ok;
            for (int i = 0; i < 1000; i++) {
                if (ratio[i] < ratio[min_index] * (1 + (double)restriction / 100.0)) {
                    ok.push_back(i);
                }
            }

            min_index = ok[rnd(ok.size())];
        }

        result.push_back(min_index);

        int del_count = 0;
        int idx = 0;
        for (int i = 0; i < uncover2.size(); i++) {
            if (count[min_index][uncover2[idx]]) {
                del_count++;
                std::swap(uncover2[idx], uncover2[uncover2.size() - del_count]);
            } else {
                idx++;
            }
        }

        uncover2.resize(uncover2.size() - del_count);
    }

    return result;
}

std::vector<int> neighbor_search(std::vector<std::vector<char>> &count, std::vector<int> &from, int priority, int restriction, int magnitude) {
    std::vector<int> result = from;

    for (int iter = 0; iter < 100; iter++) {
        int prev_size = result.size();
        std::vector<int> prev_result = result;

        int remove_num = (result.size() * magnitude + 99) / 100;
        shuffle(result);
        result.resize(result.size() - remove_num);

        std::vector<char> cover(count[0].size(), 0);
        for (int i : result) {
            for (int j = 0; j < (int)count[i].size(); j++) {
                cover[j] |= count[i][j];
            }
        }

        std::vector<int> uncover;
        for (int i = 0; i < (int)cover.size(); i++) {
            if (!cover[i]) {
                uncover.push_back(i);
            }
        }

        std::vector<int> result2 = construction(uncover, count, priority, restriction);
        result.insert(result.end(), result2.begin(), result2.end());

        if (prev_size < result.size()) {
            std::swap(result, prev_result);
        }
    }

    return result;
}

// |index| = p
// 返り値 {ans.size(),ans}
// 条件 : ans と index[i] の積集合は常に非空 (0 <= i < n)
pair<int, std::vector<int>> set_cover(int p, std::vector<std::vector<int>> &index) {
    int priority = 5;
    int restriction = 15;

    std::vector<std::vector<char>> count(1000, std::vector<char>(index.size(), 0));

    for (int i = 0; i < index.size(); i++) {
        for (int j : index[i]) {
            count[j][i] = 1;
        }
    }

    std::vector<int> ans;
    int optimal_score = p + 1;

    std::vector<int> uncover(index.size());
    for (int i = 0; i < index.size(); i++) {
        uncover[i] = i;
    }

    for (int iter = 0; iter < 100; iter++) {
        std::vector<int> x = construction(uncover, count, priority, restriction);
        int x_score = x.size();

        if (x_score < optimal_score) {
            optimal_score = x_score;
            ans = x;
        }

        if (x_score < optimal_score * (1 + (double)restriction / 100.0)) {
            x = neighbor_search(count, x, priority, restriction, 30);

            int x_score2 = x.size();
            if (x_score2 < optimal_score) {
                optimal_score = x_score2;
                ans = x;
            }
            
            std::cout << "iter " << iter << " score " << x_score << " " << x_score2 << " " << optimal_score << std::endl;
        }

        
    }

    return {ans.size(), ans};
}

/*
おきもち :
充足してる条件の数が少ないほど count[i] が少ないのも選択できるようにしたい
半分ぐらい充足してたら貪欲でいいかもしれない
*/

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    std::cin >> n;

    vvi v(n);

    rep(i, 0, n) {
        int x;
        cin >> x;
        int sz;
        cin >> sz;
        rep(j, 1, sz) {
            int val;
            cin >> val;
            v[i].emplace_back(val);
        }
    }

    /*
    rep(i, 0, n) {
        int x;
        cin >> x;
        string s;
        cin >> s;
        string s2 = "";
        rep(j, 5, s.size()) s2 += s[j];
        int sz = stoi(s2);
        cin >> s;
        s2 = "";
        rep(j, 1, s.size()) s2 += s[j];
        v[i].emplace_back(stoi(s2));
        rep(j, 1, sz) {
            int p;
            cin >> p;
            v[i].emplace_back(p);
        }
        cin >> s;
    }
    */

    auto [anssize, ans] = set_cover(n, v);
    std::cout << anssize << endl;
    for (int val : ans) std::cout << val << " ";
    std::cout << endl;
}