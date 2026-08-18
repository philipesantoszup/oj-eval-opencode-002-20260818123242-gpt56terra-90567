#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

namespace sjtu {
class int2048 {
  static const unsigned base = 10000;
  std::vector<unsigned> digit;
  bool negative;
  void trim();
  int abs_compare(const int2048 &) const;
  static int2048 abs_add(const int2048 &, const int2048 &);
  static int2048 abs_sub(const int2048 &, const int2048 &);
  static int2048 abs_multiply(const int2048 &, const int2048 &);
  static int2048 abs_divide(const int2048 &, const int2048 &);
  void multiply_unsigned(unsigned);
  void add_unsigned(unsigned);
 public:
  int2048(); int2048(long long); int2048(const std::string &); int2048(const int2048 &);
  void read(const std::string &); void print();
  int2048 &add(const int2048 &); friend int2048 add(int2048, const int2048 &);
  int2048 &minus(const int2048 &); friend int2048 minus(int2048, const int2048 &);
  int2048 operator+() const; int2048 operator-() const; int2048 &operator=(const int2048 &);
  int2048 &operator+=(const int2048 &); friend int2048 operator+(int2048, const int2048 &);
  int2048 &operator-=(const int2048 &); friend int2048 operator-(int2048, const int2048 &);
  int2048 &operator*=(const int2048 &); friend int2048 operator*(int2048, const int2048 &);
  int2048 &operator/=(const int2048 &); friend int2048 operator/(int2048, const int2048 &);
  int2048 &operator%=(const int2048 &); friend int2048 operator%(int2048, const int2048 &);
  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);
  friend bool operator==(const int2048 &, const int2048 &); friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &); friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &); friend bool operator>=(const int2048 &, const int2048 &);
};

int2048::int2048():digit(1,0),negative(false){} int2048::int2048(long long x):negative(x<0){unsigned long long y=negative?-(unsigned long long)x:(unsigned long long)x;do{digit.push_back(y%base);y/=base;}while(y);} int2048::int2048(const std::string&s){read(s);} int2048::int2048(const int2048&x)=default;
void int2048::trim(){while(digit.size()>1&&!digit.back())digit.pop_back();if(digit.size()==1&&!digit[0])negative=false;}
void int2048::read(const std::string&s){digit.clear();negative=false;std::size_t p=0;if(s[0]=='-'||s[0]=='+'){negative=s[0]=='-';p=1;}for(std::size_t e=s.size();e>p;){std::size_t b=e>=p+4?e-4:p;unsigned x=0;for(std::size_t i=b;i<e;++i)x=x*10+s[i]-'0';digit.push_back(x);e=b;}if(digit.empty())digit.push_back(0);trim();}
void int2048::print(){std::cout<<*this;} int int2048::abs_compare(const int2048&x)const{if(digit.size()!=x.digit.size())return digit.size()<x.digit.size()?-1:1;for(std::size_t i=digit.size();i-->0;)if(digit[i]!=x.digit[i])return digit[i]<x.digit[i]?-1:1;return 0;}
int2048 int2048::abs_add(const int2048&a,const int2048&b){int2048 r;r.digit.assign((a.digit.size()>b.digit.size()?a.digit.size():b.digit.size())+1,0);unsigned c=0;for(std::size_t i=0;i<r.digit.size();++i){unsigned x=c+(i<a.digit.size()?a.digit[i]:0)+(i<b.digit.size()?b.digit[i]:0);r.digit[i]=x%base;c=x/base;}r.trim();return r;}
int2048 int2048::abs_sub(const int2048&a,const int2048&b){int2048 r;r.digit=a.digit;int c=0;for(std::size_t i=0;i<r.digit.size();++i){int x=r.digit[i]-(i<b.digit.size()?b.digit[i]:0)-c;if(x<0)x+=base,c=1;else c=0;r.digit[i]=x;}r.trim();return r;}
void int2048::multiply_unsigned(unsigned x){if(!x){digit.assign(1,0);negative=false;return;}unsigned long long c=0;for(std::size_t i=0;i<digit.size();++i){unsigned long long v=(unsigned long long)digit[i]*x+c;digit[i]=v%base;c=v/base;}while(c){digit.push_back(c%base);c/=base;}}
void int2048::add_unsigned(unsigned x){for(std::size_t i=0;x;++i){if(i==digit.size())digit.push_back(0);unsigned v=digit[i]+x;digit[i]=v%base;x=v/base;}}
static void fft(std::vector<std::complex<double>>&a,bool inv){std::size_t n=a.size();for(std::size_t i=1,j=0;i<n;++i){std::size_t b=n>>1;for(;j&b;b>>=1)j^=b;j^=b;if(i<j){auto t=a[i];a[i]=a[j];a[j]=t;}}const double pi=3.14159265358979323846;for(std::size_t l=2;l<=n;l<<=1){std::complex<double>wlen(std::cos(2*pi/l),(inv?-1:1)*std::sin(2*pi/l));for(std::size_t i=0;i<n;i+=l){std::complex<double>w(1);for(std::size_t j=0;j<l/2;++j){auto u=a[i+j],v=a[i+j+l/2]*w;a[i+j]=u+v;a[i+j+l/2]=u-v;w*=wlen;}}}if(inv)for(auto&x:a)x/=(double)n;}
int2048 int2048::abs_multiply(const int2048&a,const int2048&b){if(!a.digit[0]||!b.digit[0])return int2048();if(a.digit.size()<48||b.digit.size()<48){int2048 r;r.digit.assign(a.digit.size()+b.digit.size(),0);for(std::size_t i=0;i<a.digit.size();++i){unsigned long long c=0;for(std::size_t j=0;j<b.digit.size()||c;++j){unsigned long long x=r.digit[i+j]+c+(j<b.digit.size()?(unsigned long long)a.digit[i]*b.digit[j]:0);r.digit[i+j]=x%base;c=x/base;}}r.trim();return r;}std::size_t n=1;while(n<2*(a.digit.size()+b.digit.size()))n<<=1;std::vector<std::complex<double>>x(n),y(n);for(std::size_t i=0;i<a.digit.size();++i)x[2*i]=a.digit[i]%100,x[2*i+1]=a.digit[i]/100;for(std::size_t i=0;i<b.digit.size();++i)y[2*i]=b.digit[i]%100,y[2*i+1]=b.digit[i]/100;fft(x,false);fft(y,false);for(std::size_t i=0;i<n;++i)x[i]*=y[i];fft(x,true);std::vector<unsigned>v(n);long long c=0;for(std::size_t i=0;i<n;++i){long long z=(long long)(x[i].real()+.5)+c;v[i]=z%100;c=z/100;}int2048 r;r.digit.assign((n+1)/2,0);for(std::size_t i=0;i<r.digit.size();++i)r.digit[i]=v[2*i]+(2*i+1<n?100*v[2*i+1]:0);r.trim();return r;}
int2048 int2048::abs_divide(const int2048&a,const int2048&b){if(a.abs_compare(b)<0)return int2048();unsigned norm=base/(b.digit.back()+1);int2048 x(a),y(b),r,q;x.multiply_unsigned(norm);y.multiply_unsigned(norm);q.digit.assign(x.digit.size(),0);for(std::size_t p=x.digit.size();p-->0;){r.digit.push_back(0);for(std::size_t j=r.digit.size()-1;j>0;--j)r.digit[j]=r.digit[j-1];r.digit[0]=x.digit[p];r.trim();unsigned long long hi=r.digit.back(),lo=r.digit.size()>1?r.digit[r.digit.size()-2]:0;unsigned d=(hi*base+lo)/y.digit.back();if(d>=base)d=base-1;int2048 t(y);t.multiply_unsigned(d);while(r.abs_compare(t)<0){--d;t=abs_sub(t,y);}r=abs_sub(r,t);q.digit[p]=d;}q.trim();return q;}
int2048&int2048::add(const int2048&x){return *this+=x;}int2048 add(int2048 a,const int2048&b){return a+=b;}int2048&int2048::minus(const int2048&x){return *this-=x;}int2048 minus(int2048 a,const int2048&b){return a-=b;}int2048 int2048::operator+()const{return *this;}int2048 int2048::operator-()const{int2048 r(*this);if(r.digit[0])r.negative=!r.negative;return r;}int2048&int2048::operator=(const int2048&x)=default;
int2048&int2048::operator+=(const int2048&x){if(negative==x.negative){int2048 r=abs_add(*this,x);r.negative=negative;*this=r;return *this;}int c=abs_compare(x);if(c>=0){int2048 r=abs_sub(*this,x);r.negative=negative;*this=r;}else{int2048 r=abs_sub(x,*this);r.negative=x.negative;*this=r;}return *this;}int2048 operator+(int2048 a,const int2048&b){return a+=b;}int2048&int2048::operator-=(const int2048&x){return *this+=-x;}int2048 operator-(int2048 a,const int2048&b){return a-=b;}int2048&int2048::operator*=(const int2048&x){int2048 r=abs_multiply(*this,x);r.negative=negative!=x.negative;r.trim();*this=r;return *this;}int2048 operator*(int2048 a,const int2048&b){return a*=b;}
int2048&int2048::operator/=(const int2048&x){int2048 q=abs_divide(*this,x),p=abs_multiply(q,x);bool rem=abs_compare(p)!=0;q.negative=negative!=x.negative;q.trim();if(rem&&negative!=x.negative)q-=int2048(1);*this=q;return *this;}int2048 operator/(int2048 a,const int2048&b){return a/=b;}int2048&int2048::operator%=(const int2048&x){int2048 q=*this/x;*this-=q*x;return *this;}int2048 operator%(int2048 a,const int2048&b){return a%=b;}
std::istream&operator>>(std::istream&i,int2048&x){std::string s;i>>s;x.read(s);return i;}std::ostream&operator<<(std::ostream&o,const int2048&x){if(x.negative)o<<'-';o<<x.digit.back();char b[5];for(std::size_t i=x.digit.size()-1;i-->0;){std::snprintf(b,sizeof(b),"%04u",x.digit[i]);o<<b;}return o;}bool operator==(const int2048&a,const int2048&b){return a.negative==b.negative&&a.digit==b.digit;}bool operator!=(const int2048&a,const int2048&b){return!(a==b);}bool operator<(const int2048&a,const int2048&b){if(a.negative!=b.negative)return a.negative;int c=a.abs_compare(b);return a.negative?c>0:c<0;}bool operator>(const int2048&a,const int2048&b){return b<a;}bool operator<=(const int2048&a,const int2048&b){return!(b<a);}bool operator>=(const int2048&a,const int2048&b){return!(a<b);}
}
