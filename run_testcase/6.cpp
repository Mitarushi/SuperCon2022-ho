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



tuple<uint64_t, int, vi> hamming(vi &a, vi &b) {
    uint64_t res = 0;
    int cnt = 0;
    vi rv;
    for (int i = 0; i < a.size(); i++) {
      if (a[i] != b[i]) {
        res += i;
        cnt += 1;
        if (rv.size() < 100){
          rv.pb(i);
        }
        
      }
      res *= 3141;
    }
    return {res, cnt, rv};
}

int main(){
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  int n;
  cin>>n; // 1000

  vi Ta(n),Tb(n),F(n);
  for(int i=0;i<n;i++){
    cin>>Ta[i];Ta[i]--;
  }
  for(int i=0;i<n;i++){
    cin>>Tb[i];Tb[i]--;
  }
  for(int i=0;i<n;i++){
    cin>>F[i];
  }

  int m;
  cin>>m; // 2000
  vvi accept(m,vi(n)); // accpt[i][j] := w_i が状態 j の時 受理されるか

  int ccc = 0;
  for(int i=0;i<m;i++){
    string s;
    cin>>s;
    int sz=s.size();
    for(int j=0;j<n;j++){
      int nw=j;
      for(int k=0;k<sz;k++){
        if(s[k]=='a'){
          nw=Ta[nw];
        }
        else{
          nw=Tb[nw];
        }
      }
      accept[i][j]=(F[nw]==1);
      ccc += accept[i][j];
    }
  }
  map<string,int> cnt; // これまでに同じ文字列が登場した回数

  vector<string> v(m); // accept[i][j] を並べたやつ
  vi minimum(m); // v[i] の 0 と 1 のうち数の少ない方 -> 今だけ1の個数 (L124 )
  vi minimum2(n,0); // v[i] の j 文字目のうち少ない方

  for(int i=0;i<m;i++){
    string z="";
    int one=0;
    for(int j=0;j<n;j++){
      //cout<<accept[i][j];
      z+=accept[i][j]+'0';
      one+=accept[i][j];
      minimum2[j]+=accept[i][j];
    }
    //cout << endl;
    //minimum[i]=min(one,n-one);
    minimum[i]=one;
    cnt[z]++;
    v[i]=z;
  }
  
  vector<std::tuple<int, vi ,int ,int>>res;
  std::set<std::pair<uint64_t, int>> tab;
  for (int i = 0; i < m; i++) {
    for (int j = i + 1; j < m; j++) {
      auto [h, c, vv] = hamming(accept[i], accept[j]);
      if (c != 0 && !tab.contains({h, c})) {
        res.pb({c, vv , i , j});
        tab.insert({h, c});
      }
    }
  }

  sort(all(res));

  for (int i = 0; i < 10000; i ++) {
    auto [i1,v1,i2,i3]=res[i];
    if (v1.size() >= 99) {
      break;
    }
    //cout << i << " index = " << i2<<" & "<<i3<< " : size:"<< i1 << " { ";
    //cout << i << " size:"<< i1 << " {";
    cout<< i << " " << i1 <<" ";
    for (int j = 0; j < v1.size(); j++) {
      cout << v1[j] << " ";
    }
    cout << endl;
  }
}