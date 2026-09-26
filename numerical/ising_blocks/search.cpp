// Exhaustive 2D Ising square-block influence certificate.
// Build: c++ -O3 -std=c++17 -ffp-contract=off search.cpp -o search
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cfenv>
#include <limits>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct Maximizer { double value=-1; unsigned tb=0,left=0,right=0; };

// Kronecker product of [[e^x,e^-x],[e^-x,e^x]], on each matrix row.
void transfer(std::vector<double>& z, std::vector<double>& u,
              unsigned n, unsigned ell, double ep, double em) {
  for(unsigned bit=0; bit<ell; ++bit) {
    unsigned step=1u<<bit;
    for(unsigned row=0; row<n; ++row) {
      size_t base=size_t(row)*n;
      for(unsigned start=0;start<n;start+=2*step) {
        for(unsigned k=0;k<step;++k) {
          size_t a=base+start+k,b=a+step;
          double za=z[a],zb=z[b],ua=u[a],ub=u[b];
          z[a]=ep*za+em*zb; z[b]=em*za+ep*zb;
          u[a]=ep*ua+em*ub; u[b]=em*ua+ep*ub;
        }
      }
    }
  }
}

int main(int argc,char**argv) {
  static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits==53, "IEEE binary64 is required");
  if(std::fegetround()!=FE_TONEAREST) {std::cerr<<"round-to-nearest is required\n"; return 2;}
  if(argc!=5) { std::cerr<<"usage: search weights.txt ell_min ell_max output.csv\n"; return 2; }
  std::ifstream input(argv[1]); std::string xstr; int max_power;
  input>>xstr>>max_power;
  std::vector<double> w(2*max_power+1);
  for(auto& a:w) input>>a;
  if(!input) {std::cerr<<"invalid certified weights file\n"; return 2;}
  auto weight=[&](int power){return w[power+max_power];};
  const double ep=weight(1),em=weight(-1);
  std::ofstream out(argv[4]);
  out<<"x,ell,p,c_estimate,c_lower,c_upper,top_bottom_code,left_code_minus,right_code,patterns,kappa_estimate,kappa_lower,kappa_upper\n";
  out<<std::setprecision(17);
  for(unsigned ell=std::stoul(argv[2]);ell<=std::stoul(argv[3]);++ell) {
    auto begin=std::chrono::steady_clock::now();
    if(ell<1 || ell>8 || max_power<int(ell+1)) {std::cerr<<"certificate covers 1 <= ell <= 8 with sufficient weights\n";return 2;}
    const unsigned n=1u<<ell, patterns=1u<<(2*ell-1);
    std::vector<double> init(n*n),z(n*n),u(n*n);
    std::vector<int> internal(n),nup(n);
    for(unsigned s=0;s<n;++s) {
      nup[s]=__builtin_popcount(s); internal[s]=0;
      for(unsigned r=0;r+1<ell;++r)
        internal[s]+=((s>>r)&1)==((s>>(r+1))&1)?1:-1;
      for(unsigned a=0;a<n;++a)
        init[size_t(a)*n+s]=weight(int(ell)-2*__builtin_popcount(a^s));
    }
    std::vector<Maximizer> maxima(ell);
    for(unsigned pattern=0;pattern<patterns;++pattern) {
      // Top-left exterior spin fixed to + by global spin flip.
      unsigned tb=(pattern<<1)|1u;
      unsigned top=tb&(n-1),bottom=tb>>ell;
      z=init;
      for(unsigned col=0;col<ell;++col) {
        if(col) transfer(z,u,n,ell,ep,em);
        int t=(top>>col)&1?1:-1, b=(bottom>>col)&1?1:-1;
        for(unsigned s=0;s<n;++s) {
          int first=s&1?1:-1,last=(s>>(ell-1))&1?1:-1;
          double v=weight(internal[s]+t*first+b*last);
          for(unsigned a=0;a<n;++a) {
            size_t j=size_t(a)*n+s;
            if(col) u[j]=(u[j]+double(nup[s])*z[j])*v;
            else u[j]=double(nup[s])*z[j]*v;
            z[j]*=v;
          }
        }
      }
      transfer(z,u,n,ell,ep,em);
      // Safe: partition functions are strictly positive and remain normal.
      for(size_t j=0;j<z.size();++j) u[j]/=z[j];
      for(unsigned p=0;p<ell;++p) {
        unsigned mask=1u<<p;
        for(unsigned a=0;a<n;++a) if(!(a&mask)) {
          size_t lo=size_t(a)*n,hi=size_t(a|mask)*n;
          for(unsigned b=0;b<n;++b) {
            double c=u[hi+b]-u[lo+b];
            if(c>maxima[p].value) maxima[p]={c,tb,a,b};
          }
        }
      }
    }
    // See README.md: positive-operation relative-error certificate.
    // 4096 unit-roundoff factors dominate all operation-path lengths for ell<=8.
    const double abs_error=1e-9; // rigorous bound < 1.3e-10, rounded up conservatively.
    double sum=0;for(auto m:maxima)sum+=m.value;
    double kappa=1-4*sum/(ell*ell);
    double kerr=4*abs_error/ell+1e-13;
    for(unsigned p=0;p<ell;++p) {
      auto m=maxima[p];
      out<<xstr<<','<<ell<<','<<p+1<<','<<m.value<<','<<std::nextafter(m.value-abs_error,-INFINITY)<<','<<std::nextafter(m.value+abs_error,INFINITY)<<','
         <<m.tb<<','<<m.left<<','<<m.right<<','<<patterns<<','<<kappa<<','<<std::nextafter(kappa-kerr,-INFINITY)<<','<<std::nextafter(kappa+kerr,INFINITY)<<'\n';
    }
    out.flush();
    double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    std::cerr<<"x="<<xstr<<" ell="<<ell<<" kappa="<<std::setprecision(12)<<kappa<<" c=";
    for(auto m:maxima)std::cerr<<m.value<<' ';
    std::cerr<<" time="<<elapsed<<" sec\n";
  }
}
