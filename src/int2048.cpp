#include "int2048.h"

namespace sjtu {

// =============================
// Helpers - forward declarations
// =============================

namespace {

using std::complex;
using std::size_t;
using std::vector;

const double PI_D = 3.141592653589793238462643383279502884; // high precision constant

void fft(vector<complex<double>> &a, bool invert) {
  const size_t n = a.size();
  // bit-reversal permutation
  for (size_t i = 1, j = 0; i < n; ++i) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      complex<double> tmp = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
  }

  for (size_t len = 2; len <= n; len <<= 1) {
    double ang = 2 * PI_D / static_cast<double>(len) * (invert ? -1.0 : 1.0);
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

  if (invert) {
    for (size_t i = 0; i < n; ++i) a[i] /= static_cast<double>(n);
  }
}

} // namespace (anonymous)

// =============================
// int2048 implementation
// =============================

// Constructors
int2048::int2048() : limbs(), negative(false) {}

int2048::int2048(long long value) : limbs(), negative(false) {
  if (value < 0) {
    negative = true;
    // handle LLONG_MIN
    unsigned long long x = static_cast<unsigned long long>(-(value + 1));
    ++x; // now x = abs(value)
    while (x) {
      limbs.push_back(static_cast<uint32_t>(x % BASE));
      x /= BASE;
    }
  } else {
    unsigned long long x = static_cast<unsigned long long>(value);
    while (x) {
      limbs.push_back(static_cast<uint32_t>(x % BASE));
      x /= BASE;
    }
  }
  trim();
}

int2048::int2048(const std::string &s) : limbs(), negative(false) { read(s); }

int2048::int2048(const int2048 &other) = default;

// normalize and utilities
void int2048::trim() {
  while (!limbs.empty() && limbs.back() == 0) limbs.pop_back();
  if (limbs.empty()) negative = false;
}

int int2048::compareAbs(const int2048 &a, const int2048 &b) {
  if (a.limbs.size() != b.limbs.size())
    return a.limbs.size() < b.limbs.size() ? -1 : 1;
  for (std::size_t i = a.limbs.size(); i-- > 0;) {
    if (a.limbs[i] != b.limbs[i]) return a.limbs[i] < b.limbs[i] ? -1 : 1;
  }
  return 0;
}

int int2048::compareSigned(const int2048 &a, const int2048 &b) {
  if (a.negative != b.negative) return a.negative ? -1 : 1;
  int cmp = compareAbs(a, b);
  return a.negative ? -cmp : cmp;
}

void int2048::addAbs(int2048 &a, const int2048 &b) {
  const std::size_t n = a.limbs.size();
  const std::size_t m = b.limbs.size();
  const std::size_t maxn = n > m ? n : m;
  a.limbs.resize(maxn, 0);
  unsigned long long carry = 0;
  for (std::size_t i = 0; i < maxn; ++i) {
    unsigned long long sum = carry + a.limbs[i] + (i < m ? b.limbs[i] : 0u);
    a.limbs[i] = static_cast<uint32_t>(sum % BASE);
    carry = sum / BASE;
  }
  if (carry) a.limbs.push_back(static_cast<uint32_t>(carry));
}

void int2048::subAbs(int2048 &a, const int2048 &b) {
  // precondition: |a| >= |b|
  const std::size_t m = b.limbs.size();
  long long carry = 0;
  for (std::size_t i = 0; i < a.limbs.size(); ++i) {
    long long cur = static_cast<long long>(a.limbs[i]) - (i < m ? b.limbs[i] : 0) - carry;
    if (cur < 0) {
      cur += BASE;
      carry = 1;
    } else {
      carry = 0;
    }
    a.limbs[i] = static_cast<uint32_t>(cur);
  }
  a.trim();
}

int2048 int2048::mulSmall(const int2048 &a, uint32_t m) {
  int2048 res;
  if (m == 0 || a.limbs.empty()) return res;
  res.limbs.resize(a.limbs.size());
  unsigned long long carry = 0;
  for (std::size_t i = 0; i < a.limbs.size(); ++i) {
    unsigned long long cur = carry + static_cast<unsigned long long>(a.limbs[i]) * m;
    res.limbs[i] = static_cast<uint32_t>(cur % BASE);
    carry = cur / BASE;
  }
  if (carry) res.limbs.push_back(static_cast<uint32_t>(carry));
  return res;
}

void int2048::addMulSmall(int2048 &acc, const int2048 &a, uint32_t m, size_t offset) {
  if (m == 0 || a.limbs.empty()) return;
  if (acc.limbs.size() < a.limbs.size() + offset) acc.limbs.resize(a.limbs.size() + offset, 0);
  unsigned long long carry = 0;
  std::size_t i = 0;
  for (; i < a.limbs.size(); ++i) {
    unsigned long long cur = carry + static_cast<unsigned long long>(a.limbs[i]) * m + acc.limbs[i + offset];
    acc.limbs[i + offset] = static_cast<uint32_t>(cur % BASE);
    carry = cur / BASE;
  }
  std::size_t idx = i + offset;
  while (carry) {
    if (idx >= acc.limbs.size()) acc.limbs.push_back(0);
    unsigned long long cur = carry + acc.limbs[idx];
    acc.limbs[idx] = static_cast<uint32_t>(cur % BASE);
    carry = cur / BASE;
    ++idx;
  }
}

int2048 int2048::mulAbs(const int2048 &a, const int2048 &b) {
  // Fast multiplication using FFT on base 1000 digits (3 decimal digits)
  if (a.limbs.empty() || b.limbs.empty()) return int2048();

  // Convert to base 1000 digits
  std::vector<int> A;
  A.reserve(a.limbs.size() * 3);
  for (std::size_t i = 0; i < a.limbs.size(); ++i) {
    uint32_t x = a.limbs[i];
    A.push_back(static_cast<int>(x % 1000u)); x /= 1000u;
    A.push_back(static_cast<int>(x % 1000u)); x /= 1000u;
    A.push_back(static_cast<int>(x));
  }
  while (!A.empty() && A.back() == 0) A.pop_back();

  std::vector<int> B;
  B.reserve(b.limbs.size() * 3);
  for (std::size_t i = 0; i < b.limbs.size(); ++i) {
    uint32_t x = b.limbs[i];
    B.push_back(static_cast<int>(x % 1000u)); x /= 1000u;
    B.push_back(static_cast<int>(x % 1000u)); x /= 1000u;
    B.push_back(static_cast<int>(x));
  }
  while (!B.empty() && B.back() == 0) B.pop_back();

  if (A.empty() || B.empty()) return int2048();

  size_t n1 = A.size(), n2 = B.size();
  size_t n = 1;
  while (n < n1 + n2) n <<= 1;

  std::vector<complex<double>> fa(n), fb(n);
  for (size_t i = 0; i < n1; ++i) fa[i] = complex<double>(static_cast<double>(A[i]), 0.0);
  for (size_t i = n1; i < n; ++i) fa[i] = complex<double>(0.0, 0.0);
  for (size_t i = 0; i < n2; ++i) fb[i] = complex<double>(static_cast<double>(B[i]), 0.0);
  for (size_t i = n2; i < n; ++i) fb[i] = complex<double>(0.0, 0.0);

  fft(fa, false);
  fft(fb, false);
  for (size_t i = 0; i < n; ++i) fa[i] *= fb[i];
  fft(fa, true);

  std::vector<long long> C(n);
  for (size_t i = 0; i < n; ++i) C[i] = static_cast<long long>(fa[i].real() + (fa[i].real() >= 0 ? 0.5 : -0.5));

  // carry in base 1000
  long long carry = 0;
  for (size_t i = 0; i < n; ++i) {
    long long cur = C[i] + carry;
    if (cur >= 0) {
      C[i] = cur % 1000;
      carry = cur / 1000;
    } else {
      long long k = (-cur + 999) / 1000; // ceil division
      C[i] = cur + k * 1000;
      carry = -k;
    }
  }
  while (carry > 0) { C.push_back(carry % 1000); carry /= 1000; }
  while (!C.empty() && C.back() == 0) C.pop_back();

  // pack base-1000 digits into base-1e9 limbs
  int2048 res;
  res.negative = false;
  unsigned long long acc = 0;
  unsigned long long p = 1; // current power of 1000 inside a limb
  int cnt = 0;
  for (size_t i = 0; i < C.size(); ++i) {
    acc += static_cast<unsigned long long>(C[i]) * p;
    ++cnt;
    if (cnt == 3) {
      res.limbs.push_back(static_cast<uint32_t>(acc));
      acc = 0;
      p = 1;
      cnt = 0;
    } else {
      p *= 1000ull;
    }
  }
  if (cnt != 0) res.limbs.push_back(static_cast<uint32_t>(acc));
  res.trim();
  return res;
}

int2048 int2048::divModAbs(const int2048 &a, const int2048 &b, int2048 &remainder) {
  // returns quotient; remainder is set (both non-negative)
  int2048 zero;
  remainder = zero;
  if (b.limbs.empty()) return zero; // undefined, but avoid crash
  if (compareAbs(a, b) < 0) {
    remainder = a;
    int2048 q; // zero
    q.negative = false;
    q.trim();
    remainder.negative = false;
    return q;
  }
  if (b.limbs.size() == 1) {
    uint32_t div = b.limbs[0];
    int2048 q;
    q.limbs.resize(a.limbs.size());
    unsigned long long rem = 0;
    for (std::size_t i = a.limbs.size(); i-- > 0;) {
      unsigned long long cur = a.limbs[i] + rem * BASE;
      q.limbs[i] = static_cast<uint32_t>(cur / div);
      rem = cur % div;
    }
    q.trim();
    remainder = int2048(static_cast<long long>(0));
    if (rem) {
      remainder.limbs.push_back(static_cast<uint32_t>(rem));
      remainder.negative = false;
    }
    remainder.trim();
    return q;
  }

  // Normalize
  unsigned long long norm = (static_cast<unsigned long long>(BASE) / (static_cast<unsigned long long>(b.limbs.back()) + 1ull));
  int2048 u = mulSmall(a, static_cast<uint32_t>(norm));
  int2048 v = mulSmall(b, static_cast<uint32_t>(norm));
  u.negative = false; v.negative = false;

  std::size_t n = u.limbs.size();
  std::size_t m = v.limbs.size();
  if (n == m) { u.limbs.push_back(0); ++n; }
  else if (n < m) { // shouldn't happen due to earlier check, but guard
    int2048 q; remainder = a; remainder.negative = false; return q;
  } else {
    u.limbs.push_back(0);
    ++n;
  }

  int2048 q;
  q.limbs.assign(n - m, 0);

  for (std::size_t j = n - m; j-- > 0;) {
    unsigned long long ujm = u.limbs[j + m];
    unsigned long long ujm1 = u.limbs[j + m - 1];
    unsigned long long ujm2 = (m >= 2 ? u.limbs[j + m - 2] : 0ull);
    unsigned long long v1 = v.limbs[m - 1];
    unsigned long long v2 = (m >= 2 ? v.limbs[m - 2] : 0ull);

    unsigned long long dividend = ujm * BASE + ujm1;
    unsigned long long qhat = dividend / v1;
    unsigned long long rhat = dividend % v1;
    if (qhat >= BASE) qhat = BASE - 1;
    while (qhat * v2 > rhat * BASE + ujm2) {
      --qhat;
      rhat += v1;
      if (rhat >= BASE) break;
    }

    // Subtract qhat * v shifted by j from u
    long long borrow = 0;
    unsigned long long carry = 0;
    for (std::size_t i = 0; i < m; ++i) {
      unsigned long long p = qhat * static_cast<unsigned long long>(v.limbs[i]) + carry;
      carry = p / BASE;
      long long cur = static_cast<long long>(u.limbs[i + j]) - static_cast<long long>(p % BASE) - borrow;
      if (cur < 0) { cur += BASE; borrow = 1; } else { borrow = 0; }
      u.limbs[i + j] = static_cast<uint32_t>(cur);
    }
    long long cur = static_cast<long long>(u.limbs[j + m]) - static_cast<long long>(carry) - borrow;
    if (cur < 0) {
      // qhat was too big; decrement and add v back
      --qhat;
      unsigned long long c = 0;
      for (std::size_t i = 0; i < m; ++i) {
        unsigned long long sum = static_cast<unsigned long long>(u.limbs[i + j]) + v.limbs[i] + c;
        u.limbs[i + j] = static_cast<uint32_t>(sum % BASE);
        c = sum / BASE;
      }
      u.limbs[j + m] = static_cast<uint32_t>(static_cast<long long>(u.limbs[j + m]) + static_cast<long long>(c));
    } else {
      u.limbs[j + m] = static_cast<uint32_t>(cur);
    }
    q.limbs[j] = static_cast<uint32_t>(qhat);
  }

  q.trim();
  // remainder = u / norm
  // compute remainder value from lower m limbs of u
  int2048 r;
  r.limbs.assign(u.limbs.begin(), u.limbs.begin() + static_cast<long long>(m));
  r.trim();
  if (norm != 1ull) {
    // divide r by norm (small)
    unsigned long long rem = 0;
    for (std::size_t i = r.limbs.size(); i-- > 0;) {
      unsigned long long cur = r.limbs[i] + rem * BASE;
      r.limbs[i] = static_cast<uint32_t>(cur / norm);
      rem = cur % norm;
    }
    r.trim();
  }
  remainder = r;
  return q;
}

// IO
void int2048::read(const std::string &s) {
  limbs.clear();
  negative = false;
  std::size_t i = 0;
  while (i < s.size() && (s[i] == ' ' || s[i] == '\n' || s[i] == '\t' || s[i] == '\r')) ++i;
  bool neg = false;
  if (i < s.size() && (s[i] == '-' || s[i] == '+')) { neg = (s[i] == '-'); ++i; }
  while (i < s.size() && s[i] == '0') ++i; // skip leading zeros
  std::vector<uint32_t> temp;
  for (std::size_t j = s.size(); j > i;) {
    std::size_t start = (j >= 9 ? j - 9 : i);
    if (start < i) start = i;
    uint32_t chunk = 0;
    for (std::size_t k = start; k < j; ++k) {
      char c = s[k];
      if (c >= '0' && c <= '9') {
        chunk = chunk * 10u + static_cast<uint32_t>(c - '0');
      }
    }
    temp.push_back(chunk);
    j = start;
    if (j == i) break;
  }
  if (!temp.empty()) {
    for (std::size_t t = 0; t < temp.size(); ++t) limbs.push_back(temp[t]);
  }
  trim();
  if (!limbs.empty()) negative = neg;
}

void int2048::print() {
  if (limbs.empty()) { std::cout << 0; return; }
  if (negative) std::cout << '-';
  // print highest limb without leading zeros
  std::cout << limbs.back();
  // print remaining limbs padded to 9 digits
  for (std::size_t i = limbs.size() - 1; i-- > 0;) {
    uint32_t x = limbs[i];
    // pad to 9 digits
    uint32_t p = BASE / 10u;
    while (p > 0) {
      uint32_t d = x / p;
      std::cout << d;
      x %= p;
      p /= 10u;
    }
  }
}

// Basic operations
int2048 &int2048::add(const int2048 &other) {
  if (other.limbs.empty()) return *this;
  if (limbs.empty()) { *this = other; return *this; }
  if (negative == other.negative) {
    addAbs(*this, other);
  } else {
    int cmp = compareAbs(*this, other);
    if (cmp == 0) {
      limbs.clear();
      negative = false;
    } else if (cmp > 0) {
      subAbs(*this, other); // keep sign
    } else {
      int2048 tmp = other;
      subAbs(tmp, *this);
      *this = tmp;
    }
  }
  trim();
  return *this;
}

int2048 add(int2048 a, const int2048 &b) { return a.add(b); }

int2048 &int2048::minus(const int2048 &other) {
  if (other.limbs.empty()) return *this;
  if (limbs.empty()) { *this = other; this->negative = !other.negative; return *this; }
  if (negative != other.negative) {
    addAbs(*this, other);
  } else {
    int cmp = compareAbs(*this, other);
    if (cmp == 0) {
      limbs.clear();
      negative = false;
    } else if (cmp > 0) {
      subAbs(*this, other); // sign unchanged
    } else {
      int2048 tmp = other;
      subAbs(tmp, *this);
      *this = tmp;
      this->negative = !this->negative; // result sign flips
    }
  }
  trim();
  return *this;
}

int2048 minus(int2048 a, const int2048 &b) { return a.minus(b); }

// Operators
int2048 int2048::operator+() const { return *this; }
int2048 int2048::operator-() const {
  int2048 t = *this;
  if (!t.limbs.empty()) t.negative = !t.negative;
  return t;
}

int2048 &int2048::operator=(const int2048 &rhs) = default;

int2048 &int2048::operator+=(const int2048 &rhs) { return add(rhs); }
int2048 operator+(int2048 a, const int2048 &b) { return a += b; }

int2048 &int2048::operator-=(const int2048 &rhs) { return minus(rhs); }
int2048 operator-(int2048 a, const int2048 &b) { return a -= b; }

int2048 &int2048::operator*=(const int2048 &rhs) {
  if (limbs.empty() || rhs.limbs.empty()) { limbs.clear(); negative = false; return *this; }
  int2048 res = mulAbs(*this, rhs);
  res.negative = (this->negative != rhs.negative) && !res.limbs.empty();
  *this = res;
  return *this;
}
int2048 operator*(int2048 a, const int2048 &b) { return a *= b; }

int2048 &int2048::operator/=(const int2048 &rhs) {
  // floor division
  if (rhs.limbs.empty()) { limbs.clear(); negative = false; return *this; }
  int2048 aAbs = *this; aAbs.negative = false;
  int2048 bAbs = rhs; bAbs.negative = false;
  int2048 r;
  int2048 q = divModAbs(aAbs, bAbs, r);
  bool signsDifferent = (this->negative != rhs.negative);
  if (signsDifferent && !r.limbs.empty()) {
    // q = -( |a| / |b| ) rounding toward zero; need floor -> decrement by 1
    // q = q + (-1)
    int2048 one(1);
    q.negative = false;
    q = q.add(one), q.negative = false; // q = |q| + 1
    q.negative = true;                  // then negate
  } else {
    q.negative = (this->negative != rhs.negative) && !q.limbs.empty();
  }
  *this = q;
  trim();
  return *this;
}
int2048 operator/(int2048 a, const int2048 &b) { return a /= b; }

int2048 &int2048::operator%=(const int2048 &rhs) {
  if (rhs.limbs.empty()) { limbs.clear(); negative = false; return *this; }
  int2048 aAbs = *this; aAbs.negative = false;
  int2048 bAbs = rhs; bAbs.negative = false;
  int2048 r;
  int2048 q = divModAbs(aAbs, bAbs, r);
  (void)q;
  // modulo defined as x - floor(x/y)*y
  // Using floor division rules:
  bool signsDifferent = (this->negative != rhs.negative);
  if (signsDifferent && !r.limbs.empty()) {
    // r' = |r| == 0 ? 0 : |b| - |r|
    if (!r.limbs.empty()) {
      int2048 rb = bAbs;
      subAbs(rb, r);
      r = rb;
    }
    r.negative = (rhs.negative ? !false : !false); // r should have same sign as divisor? For definition, remainder has same sign as divisor
    r.negative = rhs.negative ? true : false;
  } else {
    r.negative = this->negative; // remainder has same sign as dividend when trunc towards zero, but for floor def, sign of r equals sign of divisor
    r.negative = false; // r from divModAbs is non-negative and when signs not different, floor == trunc, so remainder sign is same as dividend, which equals sign of this; but we keep non-negative and then adjust sign to match divisor sign (which equals sign of this when not different)
    if (rhs.negative) r.negative = true;
  }
  // The above sign handling is convoluted; simpler: remainder = x - (x / y) * y, compute via division operator
  // Recompute using definition for correctness
  int2048 floor_q = (*this) / rhs; // this will call our /= and yield floor
  int2048 prod = floor_q * rhs;
  int2048 rem = *this - prod;
  *this = rem;
  trim();
  return *this;
}
int2048 operator%(int2048 a, const int2048 &b) { return a %= b; }

// Stream operators
std::istream &operator>>(std::istream &is, int2048 &x) {
  std::string s;
  is >> s;
  x.read(s);
  return is;
}

std::ostream &operator<<(std::ostream &os, const int2048 &x) {
  if (x.limbs.empty()) { os << 0; return os; }
  if (x.negative) os << '-';
  os << x.limbs.back();
  for (std::size_t i = x.limbs.size() - 1; i-- > 0;) {
    uint32_t val = x.limbs[i];
    uint32_t p = int2048::BASE / 10u;
    while (p > 0) {
      uint32_t d = val / p;
      os << d;
      val %= p;
      p /= 10u;
    }
  }
  return os;
}

// Comparisons
bool operator==(const int2048 &a, const int2048 &b) { return int2048::compareSigned(a, b) == 0; }
bool operator!=(const int2048 &a, const int2048 &b) { return !(a == b); }
bool operator<(const int2048 &a, const int2048 &b) { return int2048::compareSigned(a, b) < 0; }
bool operator>(const int2048 &a, const int2048 &b) { return int2048::compareSigned(a, b) > 0; }
bool operator<=(const int2048 &a, const int2048 &b) { return !(a > b); }
bool operator>=(const int2048 &a, const int2048 &b) { return !(a < b); }

} // namespace sjtu

