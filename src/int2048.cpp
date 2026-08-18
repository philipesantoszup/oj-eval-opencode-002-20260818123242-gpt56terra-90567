#include "int2048.h"

namespace sjtu {

int2048::int2048() : digit(1, 0), negative(false) {}

int2048::int2048(long long value) : negative(value < 0) {
  // Avoid negating LLONG_MIN in its signed type.
  unsigned long long magnitude = negative ? -(unsigned long long)value : (unsigned long long)value;
  do {
    digit.push_back((unsigned)(magnitude % base));
    magnitude /= base;
  } while (magnitude != 0);
}

int2048::int2048(const std::string &s) { read(s); }

int2048::int2048(const int2048 &other) = default;

void int2048::trim() {
  while (digit.size() > 1 && digit.back() == 0) digit.pop_back();
  if (digit.size() == 1 && digit[0] == 0) negative = false;
}

void int2048::read(const std::string &s) {
  digit.clear();
  negative = false;
  std::size_t start = 0;
  if (!s.empty() && (s[0] == '-' || s[0] == '+')) {
    negative = s[0] == '-';
    start = 1;
  }
  for (std::size_t end = s.size(); end > start;) {
    std::size_t begin = end >= start + 4 ? end - 4 : start;
    unsigned value = 0;
    for (std::size_t i = begin; i < end; ++i) value = value * 10 + (unsigned)(s[i] - '0');
    digit.push_back(value);
    end = begin;
  }
  if (digit.empty()) digit.push_back(0);
  trim();
}

void int2048::print() { std::cout << *this; }

int int2048::abs_compare(const int2048 &other) const {
  if (digit.size() != other.digit.size()) return digit.size() < other.digit.size() ? -1 : 1;
  for (std::size_t i = digit.size(); i-- > 0;)
    if (digit[i] != other.digit[i]) return digit[i] < other.digit[i] ? -1 : 1;
  return 0;
}

int2048 int2048::abs_add(const int2048 &a, const int2048 &b) {
  int2048 result;
  result.digit.assign((a.digit.size() > b.digit.size() ? a.digit.size() : b.digit.size()) + 1, 0);
  unsigned carry = 0;
  for (std::size_t i = 0; i < result.digit.size(); ++i) {
    unsigned sum = carry + (i < a.digit.size() ? a.digit[i] : 0) + (i < b.digit.size() ? b.digit[i] : 0);
    result.digit[i] = sum % base;
    carry = sum / base;
  }
  result.trim();
  return result;
}

int2048 int2048::abs_sub(const int2048 &a, const int2048 &b) {
  int2048 result;
  result.digit = a.digit;
  int borrow = 0;
  for (std::size_t i = 0; i < result.digit.size(); ++i) {
    int value = (int)result.digit[i] - (i < b.digit.size() ? (int)b.digit[i] : 0) - borrow;
    if (value < 0) value += base, borrow = 1;
    else borrow = 0;
    result.digit[i] = (unsigned)value;
  }
  result.trim();
  return result;
}

void int2048::multiply_unsigned(unsigned value) {
  if (value == 0) { digit.assign(1, 0); negative = false; return; }
  unsigned long long carry = 0;
  for (std::size_t i = 0; i < digit.size(); ++i) {
    unsigned long long current = (unsigned long long)digit[i] * value + carry;
    digit[i] = (unsigned)(current % base);
    carry = current / base;
  }
  while (carry) { digit.push_back((unsigned)(carry % base)); carry /= base; }
}

void int2048::add_unsigned(unsigned value) {
  std::size_t i = 0;
  while (value) {
    if (i == digit.size()) digit.push_back(0);
    unsigned current = digit[i] + value;
    digit[i] = current % base;
    value = current / base;
    ++i;
  }
}

static void fft(std::vector<std::complex<double> > &a, bool invert) {
  std::size_t n = a.size();
  for (std::size_t i = 1, j = 0; i < n; ++i) {
    std::size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { std::complex<double> t = a[i]; a[i] = a[j]; a[j] = t; }
  }
  const double pi = 3.1415926535897932384626433832795;
  for (std::size_t length = 2; length <= n; length <<= 1) {
    std::complex<double> root(std::cos(2 * pi / length), (invert ? -1.0 : 1.0) * std::sin(2 * pi / length));
    for (std::size_t i = 0; i < n; i += length) {
      std::complex<double> w(1);
      for (std::size_t j = 0; j < length / 2; ++j) {
        std::complex<double> u = a[i + j], v = a[i + j + length / 2] * w;
        a[i + j] = u + v; a[i + j + length / 2] = u - v; w *= root;
      }
    }
  }
  if (invert) for (std::size_t i = 0; i < n; ++i) a[i] /= (double)n;
}

int2048 int2048::abs_multiply(const int2048 &a, const int2048 &b) {
  if ((a.digit.size() == 1 && a.digit[0] == 0) || (b.digit.size() == 1 && b.digit[0] == 0)) return int2048();
  if (a.digit.size() < 48 || b.digit.size() < 48) {
    int2048 result;
    result.digit.assign(a.digit.size() + b.digit.size(), 0);
    for (std::size_t i = 0; i < a.digit.size(); ++i) {
      unsigned long long carry = 0;
      for (std::size_t j = 0; j < b.digit.size() || carry; ++j) {
        unsigned long long current = result.digit[i + j] + carry + (j < b.digit.size() ? (unsigned long long)a.digit[i] * b.digit[j] : 0);
        result.digit[i + j] = (unsigned)(current % base); carry = current / base;
      }
    }
    result.trim(); return result;
  }
  // Split each base-10000 limb into two base-100 coefficients.  This keeps
  // FFT rounding error safely below one even for the largest allowed inputs.
  std::size_t coefficient_count_a = a.digit.size() * 2, coefficient_count_b = b.digit.size() * 2;
  std::size_t n = 1;
  while (n < coefficient_count_a + coefficient_count_b) n <<= 1;
  std::vector<std::complex<double> > fa(n), fb(n);
  for (std::size_t i = 0; i < a.digit.size(); ++i) { fa[2 * i] = (double)(a.digit[i] % 100); fa[2 * i + 1] = (double)(a.digit[i] / 100); }
  for (std::size_t i = 0; i < b.digit.size(); ++i) { fb[2 * i] = (double)(b.digit[i] % 100); fb[2 * i + 1] = (double)(b.digit[i] / 100); }
  fft(fa, false); fft(fb, false);
  for (std::size_t i = 0; i < n; ++i) fa[i] *= fb[i];
  fft(fa, true);
  int2048 result;
  std::vector<unsigned> coefficient(n, 0);
  long long carry = 0;
  for (std::size_t i = 0; i < n; ++i) {
    long long current = (long long)(fa[i].real() + 0.5) + carry;
    coefficient[i] = (unsigned)(current % 100); carry = current / 100;
  }
  while (carry) { coefficient.push_back((unsigned)(carry % 100)); carry /= 100; }
  result.digit.assign((coefficient.size() + 1) / 2, 0);
  for (std::size_t i = 0; i < result.digit.size(); ++i)
    result.digit[i] = coefficient[2 * i] + (2 * i + 1 < coefficient.size() ? 100 * coefficient[2 * i + 1] : 0);
  result.trim(); return result;
}

int2048 int2048::abs_divide(const int2048 &a, const int2048 &b) {
  if (a.abs_compare(b) < 0) return int2048();
  unsigned norm = base / (b.digit.back() + 1);
  int2048 dividend(a), divisor(b), remainder, quotient;
  dividend.multiply_unsigned(norm); divisor.multiply_unsigned(norm);
  quotient.digit.assign(dividend.digit.size(), 0);
  for (std::size_t pos = dividend.digit.size(); pos-- > 0;) {
    remainder.digit.push_back(0);
    for (std::size_t j = remainder.digit.size() - 1; j > 0; --j) remainder.digit[j] = remainder.digit[j - 1];
    remainder.digit[0] = dividend.digit[pos]; remainder.trim();
    unsigned long long top = remainder.digit.back();
    unsigned long long next = remainder.digit.size() > 1 ? remainder.digit[remainder.digit.size() - 2] : 0;
    unsigned estimate = (unsigned)((top * base + next) / divisor.digit.back());
    if (estimate >= base) estimate = base - 1;
    int2048 product(divisor); product.multiply_unsigned(estimate);
    while (remainder.abs_compare(product) < 0) { --estimate; product = abs_sub(product, divisor); }
    remainder = abs_sub(remainder, product);
    quotient.digit[pos] = estimate;
  }
  quotient.trim(); return quotient;
}

int2048 &int2048::add(const int2048 &other) { return *this += other; }
int2048 add(int2048 a, const int2048 &b) { return a += b; }
int2048 &int2048::minus(const int2048 &other) { return *this -= other; }
int2048 minus(int2048 a, const int2048 &b) { return a -= b; }
int2048 int2048::operator+() const { return *this; }
int2048 int2048::operator-() const { int2048 result(*this); if (result.digit[0] != 0) result.negative = !result.negative; return result; }
int2048 &int2048::operator=(const int2048 &other) = default;

int2048 &int2048::operator+=(const int2048 &other) {
  if (negative == other.negative) { int2048 result = abs_add(*this, other); result.negative = negative; *this = result; return *this; }
  int comparison = abs_compare(other);
  if (comparison >= 0) { int2048 result = abs_sub(*this, other); result.negative = negative; *this = result; }
  else { int2048 result = abs_sub(other, *this); result.negative = other.negative; *this = result; }
  return *this;
}
int2048 operator+(int2048 a, const int2048 &b) { return a += b; }
int2048 &int2048::operator-=(const int2048 &other) { return *this += -other; }
int2048 operator-(int2048 a, const int2048 &b) { return a -= b; }
int2048 &int2048::operator*=(const int2048 &other) { int2048 result = abs_multiply(*this, other); result.negative = negative != other.negative; result.trim(); *this = result; return *this; }
int2048 operator*(int2048 a, const int2048 &b) { return a *= b; }

int2048 &int2048::operator/=(const int2048 &other) {
  int2048 quotient = abs_divide(*this, other);
  int2048 product = abs_multiply(quotient, other);
  bool has_remainder = abs_compare(product) != 0;
  quotient.negative = negative != other.negative;
  quotient.trim();
  if (has_remainder && negative != other.negative) quotient -= int2048(1);
  *this = quotient; return *this;
}
int2048 operator/(int2048 a, const int2048 &b) { return a /= b; }
int2048 &int2048::operator%=(const int2048 &other) { int2048 quotient = *this / other; *this -= quotient * other; return *this; }
int2048 operator%(int2048 a, const int2048 &b) { return a %= b; }

std::istream &operator>>(std::istream &in, int2048 &value) { std::string s; in >> s; value.read(s); return in; }
std::ostream &operator<<(std::ostream &out, const int2048 &value) {
  if (value.negative) out << '-';
  out << value.digit.back();
  char buffer[5];
  for (std::size_t i = value.digit.size() - 1; i-- > 0;) { std::snprintf(buffer, sizeof(buffer), "%04u", value.digit[i]); out << buffer; }
  return out;
}
bool operator==(const int2048 &a, const int2048 &b) { return a.negative == b.negative && a.digit == b.digit; }
bool operator!=(const int2048 &a, const int2048 &b) { return !(a == b); }
bool operator<(const int2048 &a, const int2048 &b) { if (a.negative != b.negative) return a.negative; int c = a.abs_compare(b); return a.negative ? c > 0 : c < 0; }
bool operator>(const int2048 &a, const int2048 &b) { return b < a; }
bool operator<=(const int2048 &a, const int2048 &b) { return !(b < a); }
bool operator>=(const int2048 &a, const int2048 &b) { return !(a < b); }

} // namespace sjtu
