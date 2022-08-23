//g++ 5.cpp -std=c++17 -O2 -I .
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

// |index| = p
// 返り値 {ans.size(),ans}
// 条件 : ans と index[i] の積集合は常に非空 (0 <= i < n)
pair<int,std::vector<int>> set_cover(int p,std::vector<std::vector<int>> &index){
  std::vector<int> count(1000,0); // count(sc::N_MAX,0); 残っている条件で i を選択したときいくつの条件が充足されるか
  std::vector<int> ans; // 答えを格納する集合
  std::vector<int> satisfy(p,0); // index[i] についての条件が充足されているか ans と index[i] の積集合の大きさ
  std::vector<std::vector<int>> emerge(p,vector<int>(1000,0)); // [i][j] := index[i] が j を含んでいるか
  std::vector<int> sum_index; // index[i] の和集合
  std::vector<int> count_now(1000,0); // now が含んでいるか

  /*
  焼きなまし
  現状解 := 貪欲解
  最適解 := 貪欲解

  while some_minute :
    検討解 = 現状解 + 要素いくつか
    検討解 = 検討解 - 要素いくつか (条件を保つように)

    size(検討解) と size(現状解) で焼きなまし

    if size(現状解) < size(最適解) :
      最適解 = 現状解

  return 最適解 
  */

  /*
  ビームサーチ
  vector<今の状態> current := 空 * ビーム幅
  要素を価値順にソート (充足条件の数やindex[i]の大きさなどによる？)

  while true :

  */

  /*
  {}
  |
  {1} {2} {3} ... {}
  |
  {1,2} {1,3} ... {}
  |

  */

  
  for(int i=0;i<p;i++){
    for(int val:index[i]){
      count[val]++;
      emerge[i][val]=1;
    }
  }

  for(int i=0;i<1000;i++){
    if(count[i]!=0){
      sum_index.emplace_back(i);
    }
  }

  /*
  // 貪欲解
  while(true){
    int id=-1,mx=-1;
    for(int i=0;i<1000;i++){
      if(count[i]>mx){
        mx=count[i];
        id=i;
      }
    }

    if(mx==0)break;

    ans.emplace_back(id);

    count[id]=0;

    for(int i=0;i<p;i++){
      if(emerge[i][id]==1&&satisfy[i]==0){
        satisfy[i]=1;
        for(int val:index[i]){
          count[val]--;
        }
      }
      else if(emerge[i][id]==1){
        satisfy[i]++;
      }
    }
  }
  */
  

  //ans = { 7, 14, 35, 72, 75, 87, 88, 202, 353, 406, 409, 444, 503, 560, 576, 584, 611, 620, 666, 690, 692, 728, 731, 735, 767, 794, 800, 834, 836, 837, 844, 847, 862, 912, 928 };
  //ans = { 1,2,3,4,5,6,8, 928, 406, 836, 794, 576, 35, 690, 150, 213, 837, 645, 103, 738, 182, 952, 666, 362, 786, 75, 353, 584, 125, 110, 692, 72, 87, 800, 88, 659, 834, 458, 7, 14, 170 };

  //ans = { 928, 406, 836, 794, 576, 35, 690, 150, 213, 837, 645, 103, 738, 182, 952, 666, 362, 786, 75, 353, 584, 125, 110, 692, 72, 87, 800, 88, 659, 834, 458, 7, 14, 170 };

  //ans = {46, 516, 100, 460, 75, 816, 229, 711, 339, 112, 928, 837, 622, 952, 272, 703, 35, 690, 427, 277, 213, 692, 181, 576, 360, 659, 458, 409, 202, 319, 406, 836, 560, 107, 834, 916, 190, 88};

  int ans_sz;
  cin>>ans_sz;
  ans.clear();
  rep(i,0,ans_sz){
    int id;
    cin>>id;
    ans.emplace_back(id);
  }
  
  // ans を指定する場合 ( 実験用)
  for(int id:ans){
    for(int i=0;i<p;i++){
      if(emerge[i][id]==1){
        satisfy[i]++;
      }
    }
  }


  for(int i=0;i<p;i++){
    cout<<satisfy[i]<<" ";
    if(satisfy[i]==0){
        cout<<"iiii";
        return {-1,{}};
    }
  }
  cout<<endl;
  

  /*
  // ビームサーチ
  int beam_width=250;
  std::vector<std::pair<unsigned long long,std::vector<int>>> beam_set(beam_width); // hash,beam
  std::vector<unsigned long long> zobristhash(1000);

  for(int i=0;i<1000;i++){
    zobristhash[i]=rnd(INT64_MAX);
  }

  while(true){
    std::vector<std::tuple<int,unsigned long long,std::vector<int>>> tmp; // スコア,hash,集合
    std::set<std::pair<int,unsigned long long>> seen; // score_bb,hash_bb

    cout<<beam_set.size()<<endl;

    for(auto &[hash,bb]:beam_set){
      std::vector<int> is_it_satisfy(p,0);
      for(int val:bb){
        for(int j=0;j<p;j++){
          if(emerge[j][val]==1){
            is_it_satisfy[j]++;
          }
        }
      }

      for(int i=0;i<1000;i++){
        bb.emplace_back(i);

        for(int j=0;j<p;j++){
          if(emerge[j][i]==1){
            is_it_satisfy[j]++;
          }
        }
        
        int score_bb=0; // 後でやtta
        for(int j=0;j<p;j++){
          if(is_it_satisfy[j]!=0)score_bb++;
        }

        for(int j=0;j<p;j++){
          if(emerge[j][i]==1){
            is_it_satisfy[j]--;
          }
        }

        score_bb=1000-score_bb;

        unsigned long long hash_bb=(hash+zobristhash[i]); // 後でやtta

        if(seen.count({score_bb,hash_bb})!=0){
          bb.pop_back();
          continue;
        }

        seen.insert({score_bb,hash_bb});
        tmp.emplace_back(score_bb,hash_bb,bb);
        bb.pop_back();
      }
    }

    sort(tmp.begin(),tmp.end());
    beam_set.clear();

    cout<<get<2>(tmp[0]).size()<<" : "<<get<0>(tmp[0])<<endl;

    bool finished=0;
    while(beam_set.size()<beam_width){
      for(auto [score_bb,hash_bb,bb]:tmp){
        //cout<<bb.size()<<" ";
          beam_set.emplace_back(hash_bb,bb);
          if(score_bb==1000-p){
            ans=bb;
            finished=1;
            break;
          }
          if(beam_set.size()==beam_width)break;
      }
    }

    //cout<<endl;

    if(finished==1){
      break;
    }
  }
  */

  /*
  for(int i=0;i<p;i++){
    std::cout<<satisfy[i]<<" ";
  }
  std::cout<<endl;

  cout<<"initial_state = "<<ans.size()<<endl;
  for(int val:ans){
    cout<<val<<" ";
  }
  cout<<endl;
  
  //return {ans.size(),ans};

  // ans は貪欲解

  std::vector<int> now=ans;
  for(int val:now)count_now[val]++;

  double start_temp=21,end_temp=1;
  int sz=sum_index.size();
  int lp=10000;

  cout<<sz<<endl;

  int ch=0;

  for(int loop=0;loop<lp;loop++){
    double temp = start_temp + (end_temp - start_temp) * (double)loop/1000;
    int delta=1+(20*(lp-loop))/lp;

    std::vector<int> nxt=now;
    std::vector<int> satisfy_nxt=satisfy;
    std::vector<int> add;

    for(int i=0;i<delta;i++){
      for(int j=0;j<100;j++){
        int rd=rnd(sz);
        int r=sum_index[rd];
        if(count_now[r]!=0)continue;
        nxt.emplace_back(r);
        add.emplace_back(r);

        for(int k=0;k<p;k++){
          if(emerge[k][r]==1){
            satisfy_nxt[k]++;
          }
        }
        break;
      }
    }

    for(int i=0;i<500;i++){
      int sz2=nxt.size();
      int rd1=rnd(sz2);
      swap(nxt[rd1],nxt[sz2-1]);

      // 一番後ろを消してもいいか
      int consider=nxt[sz2-1];

      bool flg=1;
      for(int j=0;j<p;j++){
        if(emerge[j][consider]==1){
          if(satisfy_nxt[j]==1){
            flg=0;break;
          }
        }
      }

      if(flg==0)continue;

      nxt.pop_back();
      for(int j=0;j<p;j++){
        if(emerge[j][consider]==1){
          satisfy_nxt[j]--;
        }
      }
    }


    double prob=100;
    if(now.size()<nxt.size())prob=exp((now.size()-nxt.size())/temp);

    if(prob > (rnd(1000))/(double)1000){
      now = nxt;
      satisfy = satisfy_nxt;

      ch++;

      //cout<<now.size()<<endl;

      if(now.size()<ans.size()){
        ans=now;
        cout<<ans.size()<<endl;
      }
    }
  }

  cout<<"now = "<<now.size()<<" : ch = "<<ch<<endl;
  */

  return {ans.size(),ans};
}

/*
おきもち : 
充足してる条件の数が少ないほど count[i] が少ないのも選択できるようにしたい
半分ぐらい充足してたら貪欲でいいかもしれない
*/

int main(){
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  int n;
  cin>>n;
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
  rep(i,0,n){
    int x;
    cin>>x;
    string s;
    cin>>s;
    string s2="";
    rep(j,5,s.size())s2+=s[j];
    int sz=stoi(s2);
    cin>>s;
    s2="";
    rep(j,1,s.size())s2+=s[j];
    v[i].emplace_back(stoi(s2));
    rep(j,1,sz){
      int p;
      cin>>p;
      v[i].emplace_back(p);
    }
    cin>>s;
  }
  */

  auto [anssize,ans]=set_cover(n,v);
  std::cout<<anssize<<endl;
  for(int val:ans)std::cout<<val<<" ";
  std::cout<<endl;
}