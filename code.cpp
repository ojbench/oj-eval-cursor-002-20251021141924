#include <bits/stdc++.h>
using namespace std;

// Merge of header and implementation for OJ submission
namespace sjtu {
class int2048 {
private:
  static const uint32_t BASE = 1000000000u;
  static const uint32_t BASE_DIGS = 9u;
  std::vector<uint32_t> limbs;
  bool negative = false;

  void trim();
  static int compareAbs(const int2048 &a, const int2048 &b);
  static int compareSigned(const int2048 &a, const int2048 &b);
  static void addAbs(int2048 &a, const int2048 &b);
  static void subAbs(int2048 &a, const int2048 &b);
  static int2048 mulAbs(const int2048 &a, const int2048 &b);
  static int2048 mulSmall(const int2048 &a, uint32_t m);
  static void addMulSmall(int2048 &acc, const int2048 &a, uint32_t m, size_t offset);
  static int2048 divModAbs(const int2048 &a, const int2048 &b, int2048 &remainder);

public:
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);

  void read(const std::string &);
  void print();

  int2048 &add(const int2048 &);
  friend int2048 add(int2048, const int2048 &);

  int2048 &minus(const int2048 &);
  friend int2048 minus(int2048, const int2048 &);

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);

  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);

  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);

  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);

  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);

  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);

  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);

  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};

} // namespace sjtu

namespace {
using std::complex; using std::size_t; using std::vector; const double PI_D = 3.141592653589793238462643383279502884;
void fft(vector<complex<double>> &a, bool invert) {
  const size_t n = a.size();
  for (size_t i = 1, j = 0; i < n; ++i) { size_t bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; if (i < j) { auto t=a[i]; a[i]=a[j]; a[j]=t; } }
  for (size_t len = 2; len <= n; len <<= 1) {
    double ang = 2 * PI_D / (double)len * (invert ? -1.0 : 1.0);
    complex<double> wlen = std::polar(1.0, ang);
    for (size_t i = 0; i < n; i += len) {
      complex<double> w(1.0, 0.0);
      for (size_t j = 0; j < (len >> 1); ++j) {
        complex<double> u = a[i + j];
        complex<double> v = a[i + j + (len >> 1)] * w;
        a[i + j] = u + v;
        a[i + j + (len >> 1)] = u - v;
        w *= wlen;
      }
    }
  }
  if (invert) for (size_t i = 0; i < n; ++i) a[i] /= (double)n;
}
}

namespace sjtu {

int2048::int2048() : limbs(), negative(false) {}
int2048::int2048(long long value) : limbs(), negative(false) {
  if (value < 0) { negative = true; unsigned long long x = (unsigned long long)(-(value + 1)); ++x; while (x) { limbs.push_back((uint32_t)(x % BASE)); x /= BASE; } }
  else { unsigned long long x = (unsigned long long)value; while (x) { limbs.push_back((uint32_t)(x % BASE)); x /= BASE; } }
  trim();
}
int2048::int2048(const std::string &s) : limbs(), negative(false) { read(s); }
int2048::int2048(const int2048 &other) = default;

void int2048::trim() { while (!limbs.empty() && limbs.back() == 0) limbs.pop_back(); if (limbs.empty()) negative = false; }
int int2048::compareAbs(const int2048 &a, const int2048 &b) {
  if (a.limbs.size() != b.limbs.size()) return a.limbs.size() < b.limbs.size() ? -1 : 1;
  for (size_t i = a.limbs.size(); i-- > 0;) { if (a.limbs[i] != b.limbs[i]) return a.limbs[i] < b.limbs[i] ? -1 : 1; }
  return 0;
}
int int2048::compareSigned(const int2048 &a, const int2048 &b) {
  if (a.negative != b.negative) return a.negative ? -1 : 1; int cmp = compareAbs(a,b); return a.negative ? -cmp : cmp;
}
void int2048::addAbs(int2048 &a, const int2048 &b) {
  size_t n=a.limbs.size(), m=b.limbs.size(), maxn=n>m?n:m; a.limbs.resize(maxn,0); unsigned long long carry=0;
  for (size_t i=0;i<maxn;++i){ unsigned long long sum=carry+a.limbs[i]+(i<m?b.limbs[i]:0u); a.limbs[i]=(uint32_t)(sum%BASE); carry=sum/BASE; }
  if (carry) a.limbs.push_back((uint32_t)carry);
}
void int2048::subAbs(int2048 &a, const int2048 &b) {
  size_t m=b.limbs.size(); long long carry=0; for (size_t i=0;i<a.limbs.size();++i){ long long cur=(long long)a.limbs[i] - (i<m?b.limbs[i]:0) - carry; if (cur<0){ cur+=BASE; carry=1; } else carry=0; a.limbs[i]=(uint32_t)cur; } a.trim();
}
int2048 int2048::mulSmall(const int2048 &a, uint32_t m){ int2048 res; if (m==0||a.limbs.empty()) return res; res.limbs.resize(a.limbs.size()); unsigned long long carry=0; for(size_t i=0;i<a.limbs.size();++i){ unsigned long long cur=carry+(unsigned long long)a.limbs[i]*m; res.limbs[i]=(uint32_t)(cur%BASE); carry=cur/BASE; } if(carry) res.limbs.push_back((uint32_t)carry); return res; }
void int2048::addMulSmall(int2048 &acc, const int2048 &a, uint32_t m, size_t offset){ if(m==0||a.limbs.empty()) return; if(acc.limbs.size()<a.limbs.size()+offset) acc.limbs.resize(a.limbs.size()+offset,0); unsigned long long carry=0; size_t i=0; for(;i<a.limbs.size();++i){ unsigned long long cur=carry+(unsigned long long)a.limbs[i]*m+acc.limbs[i+offset]; acc.limbs[i+offset]=(uint32_t)(cur%BASE); carry=cur/BASE; } size_t idx=i+offset; while(carry){ if(idx>=acc.limbs.size()) acc.limbs.push_back(0); unsigned long long cur=carry+acc.limbs[idx]; acc.limbs[idx]=(uint32_t)(cur%BASE); carry=cur/BASE; ++idx; } }

int2048 int2048::mulAbs(const int2048 &a, const int2048 &b){ if(a.limbs.empty()||b.limbs.empty()) return int2048(); vector<int>A; A.reserve(a.limbs.size()*3); for(size_t i=0;i<a.limbs.size();++i){ uint32_t x=a.limbs[i]; A.push_back((int)(x%1000u)); x/=1000u; A.push_back((int)(x%1000u)); x/=1000u; A.push_back((int)(x)); } while(!A.empty()&&A.back()==0) A.pop_back(); vector<int>B; B.reserve(b.limbs.size()*3); for(size_t i=0;i<b.limbs.size();++i){ uint32_t x=b.limbs[i]; B.push_back((int)(x%1000u)); x/=1000u; B.push_back((int)(x%1000u)); x/=1000u; B.push_back((int)(x)); } while(!B.empty()&&B.back()==0) B.pop_back(); if(A.empty()||B.empty()) return int2048(); size_t n1=A.size(), n2=B.size(), n=1; while(n<n1+n2) n<<=1; vector<complex<double>> fa(n), fb(n); for(size_t i=0;i<n1;++i) fa[i]=complex<double>((double)A[i],0.0); for(size_t i=n1;i<n;++i) fa[i]=complex<double>(0.0,0.0); for(size_t i=0;i<n2;++i) fb[i]=complex<double>((double)B[i],0.0); for(size_t i=n2;i<n;++i) fb[i]=complex<double>(0.0,0.0); fft(fa,false); fft(fb,false); for(size_t i=0;i<n;++i) fa[i]*=fb[i]; fft(fa,true); vector<long long>C(n); for(size_t i=0;i<n;++i) C[i]=(long long)(fa[i].real() + (fa[i].real()>=0?0.5:-0.5)); long long carry=0; for(size_t i=0;i<n;++i){ long long cur=C[i]+carry; if(cur>=0){ C[i]=cur%1000; carry=cur/1000; } else { long long k=(-cur+999)/1000; C[i]=cur+k*1000; carry=-k; } } while(carry>0){ C.push_back(carry%1000); carry/=1000; } while(!C.empty()&&C.back()==0) C.pop_back(); int2048 res; res.negative=false; unsigned long long acc=0, p=1; int cnt=0; for(size_t i=0;i<C.size();++i){ acc += (unsigned long long)C[i]*p; ++cnt; if(cnt==3){ res.limbs.push_back((uint32_t)acc); acc=0; p=1; cnt=0; } else { p*=1000ull; } } if(cnt!=0) res.limbs.push_back((uint32_t)acc); res.trim(); return res; }

int2048 int2048::divModAbs(const int2048 &a, const int2048 &b, int2048 &remainder){ int2048 zero; remainder=zero; if(b.limbs.empty()) return zero; if(compareAbs(a,b)<0){ remainder=a; int2048 q; q.negative=false; q.trim(); remainder.negative=false; return q; } if(b.limbs.size()==1){ uint32_t div=b.limbs[0]; int2048 q; q.limbs.resize(a.limbs.size()); unsigned long long rem=0; for(size_t i=a.limbs.size(); i-- > 0;){ unsigned long long cur=a.limbs[i] + rem*BASE; q.limbs[i]=(uint32_t)(cur/div); rem=cur%div; } q.trim(); remainder=int2048((long long)0); if(rem){ remainder.limbs.push_back((uint32_t)rem); remainder.negative=false; } remainder.trim(); return q; }
  unsigned long long norm = ((unsigned long long)BASE) / ((unsigned long long)b.limbs.back() + 1ull);
  int2048 u = mulSmall(a, (uint32_t)norm); int2048 v = mulSmall(b, (uint32_t)norm); u.negative=false; v.negative=false; size_t n=u.limbs.size(), m=v.limbs.size(); if(n==m){ u.limbs.push_back(0); ++n; } else if(n<m){ int2048 q; remainder=a; remainder.negative=false; return q; } else { u.limbs.push_back(0); ++n; }
  int2048 q; q.limbs.assign(n-m,0);
  for(size_t j=n-m; j-- > 0;){ unsigned long long ujm=u.limbs[j+m]; unsigned long long ujm1=u.limbs[j+m-1]; unsigned long long ujm2=(m>=2?u.limbs[j+m-2]:0ull); unsigned long long v1=v.limbs[m-1]; unsigned long long v2=(m>=2?v.limbs[m-2]:0ull); unsigned long long dividend=ujm*BASE + ujm1; unsigned long long qhat=dividend / v1; unsigned long long rhat=dividend % v1; if(qhat>=BASE) qhat=BASE-1; while(qhat*v2 > rhat*BASE + ujm2){ --qhat; rhat += v1; if(rhat>=BASE) break; }
    long long borrow=0; unsigned long long carry=0; for(size_t i=0;i<m;++i){ unsigned long long p = qhat*(unsigned long long)v.limbs[i] + carry; carry = p / BASE; long long cur = (long long)u.limbs[i+j] - (long long)(p % BASE) - borrow; if(cur<0){ cur+=BASE; borrow=1; } else borrow=0; u.limbs[i+j]=(uint32_t)cur; }
    long long cur = (long long)u.limbs[j+m] - (long long)carry - borrow; if(cur<0){ --qhat; unsigned long long c=0; for(size_t i=0;i<m;++i){ unsigned long long sum=(unsigned long long)u.limbs[i+j] + v.limbs[i] + c; u.limbs[i+j]=(uint32_t)(sum%BASE); c = sum/BASE; } u.limbs[j+m] = (uint32_t)((long long)u.limbs[j+m] + (long long)c); } else { u.limbs[j+m]=(uint32_t)cur; }
    q.limbs[j]=(uint32_t)qhat; }
  q.trim(); int2048 r; r.limbs.assign(u.limbs.begin(), u.limbs.begin() + (long long)m); r.trim(); if(norm!=1ull){ unsigned long long rem=0; for(size_t i=r.limbs.size(); i-- > 0;){ unsigned long long cur=r.limbs[i] + rem*BASE; r.limbs[i]=(uint32_t)(cur / norm); rem = cur % norm; } r.trim(); }
  remainder=r; return q; }

void int2048::read(const std::string &s){ limbs.clear(); negative=false; size_t i=0; while(i<s.size() && (s[i]==' '||s[i]=='\n'||s[i]=='\t'||s[i]=='\r')) ++i; bool neg=false; if(i<s.size() && (s[i]=='-'||s[i]=='+')){ neg=(s[i]=='-'); ++i; } while(i<s.size() && s[i]=='0') ++i; vector<uint32_t> temp; for(size_t j=s.size(); j>i;){ size_t start=(j>=9?j-9:i); if(start<i) start=i; uint32_t chunk=0; for(size_t k=start;k<j;++k){ char c=s[k]; if(c>='0'&&c<='9'){ chunk = chunk*10u + (uint32_t)(c-'0'); } } temp.push_back(chunk); j=start; if(j==i) break; } if(!temp.empty()) for(size_t t=0;t<temp.size();++t) limbs.push_back(temp[t]); trim(); if(!limbs.empty()) negative=neg; }
void int2048::print(){ if(limbs.empty()){ std::cout<<0; return; } if(negative) std::cout<<'-'; std::cout<<limbs.back(); for(size_t i=limbs.size()-1; i-- > 0;){ uint32_t x=limbs[i]; uint32_t p=BASE/10u; while(p>0){ uint32_t d=x/p; std::cout<<d; x%=p; p/=10u; } } }

int2048 &int2048::add(const int2048 &other){ if(other.limbs.empty()) return *this; if(limbs.empty()){ *this=other; return *this; } if(negative==other.negative){ addAbs(*this, other); } else { int cmp=compareAbs(*this, other); if(cmp==0){ limbs.clear(); negative=false; } else if(cmp>0){ subAbs(*this, other); } else { int2048 tmp=other; subAbs(tmp, *this); *this=tmp; } } trim(); return *this; }
int2048 add(int2048 a, const int2048 &b){ return a.add(b); }
int2048 &int2048::minus(const int2048 &other){ if(other.limbs.empty()) return *this; if(limbs.empty()){ *this=other; this->negative=!other.negative; return *this; } if(negative!=other.negative){ addAbs(*this, other); } else { int cmp=compareAbs(*this, other); if(cmp==0){ limbs.clear(); negative=false; } else if(cmp>0){ subAbs(*this, other); } else { int2048 tmp=other; subAbs(tmp, *this); *this=tmp; this->negative=!this->negative; } } trim(); return *this; }
int2048 minus(int2048 a, const int2048 &b){ return a.minus(b); }

int2048 int2048::operator+() const { return *this; }
int2048 int2048::operator-() const { int2048 t=*this; if(!t.limbs.empty()) t.negative=!t.negative; return t; }
int2048 &int2048::operator=(const int2048 &rhs) = default;
int2048 &int2048::operator+=(const int2048 &rhs){ return add(rhs); }
int2048 operator+(int2048 a, const int2048 &b){ return a+=b; }
int2048 &int2048::operator-=(const int2048 &rhs){ return minus(rhs); }
int2048 operator-(int2048 a, const int2048 &b){ return a-=b; }
int2048 &int2048::operator*=(const int2048 &rhs){ if(limbs.empty()||rhs.limbs.empty()){ limbs.clear(); negative=false; return *this; } int2048 res=mulAbs(*this, rhs); res.negative=(this->negative!=rhs.negative)&&!res.limbs.empty(); *this=res; return *this; }
int2048 operator*(int2048 a, const int2048 &b){ return a*=b; }
int2048 &int2048::operator/=(const int2048 &rhs){ if(rhs.limbs.empty()){ limbs.clear(); negative=false; return *this; } int2048 aAbs=*this; aAbs.negative=false; int2048 bAbs=rhs; bAbs.negative=false; int2048 r; int2048 q=divModAbs(aAbs,bAbs,r); bool signsDifferent=(this->negative!=rhs.negative); if(signsDifferent && !r.limbs.empty()){ int2048 one(1); q.negative=false; q = q.add(one); q.negative=true; } else { q.negative=(this->negative!=rhs.negative)&&!q.limbs.empty(); } *this=q; trim(); return *this; }
int2048 operator/(int2048 a, const int2048 &b){ return a/=b; }
int2048 &int2048::operator%=(const int2048 &rhs){ int2048 floor_q = (*this) / rhs; int2048 prod = floor_q * rhs; int2048 rem = *this - prod; *this = rem; trim(); return *this; }
int2048 operator%(int2048 a, const int2048 &b){ return a%=b; }

std::istream &operator>>(std::istream &is, int2048 &x){ std::string s; is>>s; x.read(s); return is; }
std::ostream &operator<<(std::ostream &os, const int2048 &x){ if(x.limbs.empty()){ os<<0; return os; } if(x.negative) os<<'-'; os<<x.limbs.back(); for(size_t i=x.limbs.size()-1; i-- > 0;){ uint32_t val=x.limbs[i]; uint32_t p=int2048::BASE/10u; while(p>0){ uint32_t d=val/p; os<<d; val%=p; p/=10u; } } return os; }

bool operator==(const int2048 &a, const int2048 &b){ return int2048::compareSigned(a,b)==0; }
bool operator!=(const int2048 &a, const int2048 &b){ return !(a==b); }
bool operator<(const int2048 &a, const int2048 &b){ return int2048::compareSigned(a,b)<0; }
bool operator>(const int2048 &a, const int2048 &b){ return int2048::compareSigned(a,b)>0; }
bool operator<=(const int2048 &a, const int2048 &b){ return !(a>b); }
bool operator>=(const int2048 &a, const int2048 &b){ return !(a<b); }

} // namespace sjtu

int main(){
  using sjtu::int2048;
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);
  int n; if(!(cin>>n)) return 0; while(n--){ string op; cin>>op; if(op=="print"){ int2048 x; cin>>x; cout<<x<<"\n"; } else if(op=="add"){ int2048 a,b; cin>>a>>b; cout<<(a+b)<<"\n"; } else if(op=="sub"){ int2048 a,b; cin>>a>>b; cout<<(a-b)<<"\n"; } else if(op=="mul"){ int2048 a,b; cin>>a>>b; cout<<(a*b)<<"\n"; } else if(op=="div"){ int2048 a,b; cin>>a>>b; cout<<(a/b)<<"\n"; } else if(op=="mod"){ int2048 a,b; cin>>a>>b; cout<<(a%b)<<"\n"; } }
  return 0;
}
