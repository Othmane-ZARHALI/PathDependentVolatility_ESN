// Coupled 64/128/512-step stock checks on identical fine OU/stress/equity paths.
#include <cmath>
#include <vector>
#include <thread>
#include <cstdint>
static inline uint64_t hash64(uint64_t x){x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);return x^(x>>31);}
static inline double gaussian(uint64_t seed,uint64_t i){uint64_t a=hash64(seed+2*i*UINT64_C(0x9e3779b97f4a7c15)),b=hash64(seed+(2*i+1)*UINT64_C(0x9e3779b97f4a7c15));double u=(double(a>>11)+.5)*0x1.0p-53,v=(double(b>>11)+.5)*0x1.0p-53;return std::sqrt(-2*std::log(u))*std::cos(6.2831853071795864769*v);}
static inline double softplus(double z){return z>0?z+std::log1p(std::exp(-z)):std::log1p(std::exp(z));}
extern "C" void coupled_stock(int nt,int np,const double*features,const double*dw,
 const double*time,const double*dt,const int*seg,const int*mask64,const int*mask128,
 const double*lognorm,const double*readout,const double*x,uint64_t seed,double*out){
 double L=std::log(365.),lhf=L*x[21],lhs=lhf+(L-lhf)*x[18],aS=std::log(2.)*365/std::exp(lhs),aF=std::log(2.)*365/std::exp(lhf),bS=5*x[16],bF=5*x[20],w=x[17],r=x[22],floor=std::pow(.05+.07*x[19],2);
 std::vector<double>dS(nt),dF(nt),nS(nt),nF(nt),qS(nt),qF(nt);
 for(int t=0;t<nt;++t){dS[t]=std::exp(-aS*dt[t]);dF[t]=std::exp(-aF*dt[t]);nS[t]=std::sqrt(-std::expm1(-2*aS*dt[t]));nF[t]=std::sqrt(-std::expm1(-2*aF*dt[t]));qS[t]=-std::expm1(-2*aS*time[t]);qF[t]=-std::expm1(-2*aF*time[t]);}
 auto worker=[&](int begin,int end){int n=end-begin;std::vector<double>zS(n,0),zF(n,0),logs(n*3,0),v64(n,0),v128(n,0);
  for(int t=0;t<nt;++t)for(int i=begin;i<end;++i){int j=i-begin;uint64_t k=uint64_t(t)*np+i;const double*f=features+3*k;double eta=readout[0]*f[0]+readout[1]*std::pow(f[1]-readout[2],2)+readout[3]*std::pow(f[2]-readout[4],2);double lg=2*(readout[5]-softplus(readout[6]*(readout[5]-eta))/readout[6]);double old=std::exp(lg-lognorm[t]),A=std::exp(bS*zS[j]-.5*bS*bS*qS[t]),B=std::exp(bF*zF[j]-.5*bF*bF*qF[t]),Q=(1-w)*old+w*((1-r)*A+r*B),v=floor+(.2-floor)*x[seg[t]]*Q;
   if(mask64[t])v64[j]=v;if(mask128[t])v128[j]=v;
   double vv[3]={v64[j],v128[j],v};for(int l=0;l<3;++l)logs[j*3+l]+=-.5*vv[l]*dt[t]+std::sqrt(vv[l])*dw[k];
   zS[j]=dS[t]*zS[j]+nS[t]*gaussian(seed,2*k);zF[j]=dF[t]*zF[j]+nF[t]*gaussian(seed,2*k+1);
  }
  for(int i=begin;i<end;++i)for(int l=0;l<3;++l)out[l*np+i]=std::exp(logs[(i-begin)*3+l]);
 };
 std::vector<std::thread> threads;for(int j=0;j<4;++j)threads.emplace_back(worker,np*j/4,np*(j+1)/4);for(auto&t:threads)t.join();
}
