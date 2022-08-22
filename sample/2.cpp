// set PATH=C:\msys64\mingw64\bin;C:\msys64\usr\local\bin;C:\msys64\usr\bin;C:\msys64\bin;%PATH%
// set MSYSTEM=MINGW64
// bash
//g++ 2.cpp -std=c++14 -O2 -I .
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
    }
  }

  map<string,int> cnt; // これまでに同じ文字列が登場した回数

  vector<string> v(m); // accept[i][j] を並べたやつ
  vi minimum(m); // v[i] の 0 と 1 のうち数の少ない方
  vi minimum2(n,0); // v[i] の j 文字目のうち少ない方

  for(int i=0;i<m;i++){
    string z="";
    int one=0;
    for(int j=0;j<n;j++){
      cout<<accept[i][j];
      z+=accept[i][j]+'0';
      one+=accept[i][j];
      minimum2[j]+=accept[i][j];
    }
    cout<<" : "<<cnt[z]<<" :: "<<min(one,n-one)<<endl;
    minimum[i]=min(one,n-one);
    cnt[z]++;
    v[i]=z;
  }

  for(int i=0;i<n;i++){
    minimum2[i]=min(minimum2[i],m-minimum2[i]);
  }

  cout<<"compile 17:44"<<endl;
  sort(minimum2.begin(),minimum2.end());
  for(int val:minimum2)cout<<val<<" ";
  cout<<endl;

  sort(minimum.begin(),minimum.end());
  for(int val:minimum)cout<<val<<" ";
  cout<<endl;

  int run;
  //cin>>run;

  run=5;

  int g=100; // | Q' |

  // 01 両方拾えない文字列の期待値
  vector<long double> failprob(n+1,0);
  for(int k=1;k<n;k++){
    int zero=k,one=n-k;
    // zero C r + one C r
    // /
    // 1000 C r

    ld x;
    ld probzero=0,probone=0;

    if(zero>=g){
      x=0;
      rep(i,0,g){
        // *= zero-i
        x+=log(zero-i);
        // /= (i+1)
        x-=log(1000-i);
      }
      x=exp(x);
      probzero=x;
    }

    if(one>=g){
      x=0;
      rep(i,0,g){
        // *= zero-i
        x+=log(one-i);
        // /= (i+1)
        x-=log(1000-i);
      }
      x=exp(x);
      probone=x;
    }

    // 0 を k 個含むやつから全部 0 or 全部 1 を採用する確率
    failprob[k]=probzero+probone;
  }

  ld fail_expected=0;
  for(int i=0;i<m;i++){
    int zero=minimum[i];
    fail_expected+=failprob[zero];
  }

  // g 個ランダムに選んで 条件をまだ満たしていない文字列の期待値
  cout<<fail_expected<<endl;

  for(int t=0;t<run;t++){
    vector<int> cnt(n,0);
    vector<int> use;
    while(use.size()<g){
      int id=range_rnd(0,n-1);
      if(cnt[id]==0){
        cnt[id]=1;
        use.emplace_back(id);
      }
    }

    vector<string> v2(m,"");
    for(int i=0;i<m;i++){
      for(int j=0;j<g;j++){
        v2[i]+=accept[i][use[j]]+'0';
      }
    }

    int fail_pair=0; // 失敗ペア
    for(int i=0;i<m;i++){
      for(int j=i+1;j<m;j++){
        if(v[i]!=v[j]&&v2[i]==v2[j]&&min(minimum[i],minimum[j])>=50){
          cout<<i<<" "<<j<<" "<<v2[i]<<" "<<v2[j]<<endl;
          fail_pair++;
        }
      }
    }

    cout<<fail_pair<<endl;
  }
}