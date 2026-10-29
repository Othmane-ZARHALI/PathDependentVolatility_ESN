// Dedicated independent OU sensitivity state. Historical cap code stays unchanged.
#include "cap_tail_kernel.cpp"
#include <cstdint>

static inline uint64_t mix64(uint64_t x){
 x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);
 x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);return x^(x>>31);
}
static inline double innovation(uint64_t seed,uint64_t index){
 uint64_t a=mix64(seed+2*index*UINT64_C(0x9e3779b97f4a7c15));
 uint64_t b=mix64(seed+(2*index+1)*UINT64_C(0x9e3779b97f4a7c15));
 double u=(double(a>>11)+.5)*0x1.0p-53,v=(double(b>>11)+.5)*0x1.0p-53;
 return std::sqrt(-2*std::log(u))*std::cos(6.2831853071795864769*v);
}
extern "C" void innovations(uint64_t seed,int n,double*out){for(int i=0;i<n;++i)out[i]=innovation(seed,i);}

// Coefficients are a Chebyshev expansion in A on [0.5,2], row per old state.
// This integrates over the conditional Gaussian Z, not over its terminal draw.
extern "C" void integrate_response(int n,int nc,int nq,const double*c,const double*mu,
 double sd,const double*z,const double*w,double sensitivity,double*out){
 for(int i=0;i<n;++i){double value=0;
  for(int k=0;k<nq;++k){double A=std::exp(sensitivity*std::tanh(mu[i]+sd*z[k]));
   double x=(2*A-2.5)/1.5,b1=0,b2=0;const double*cc=c+i*nc;
   for(int j=nc-1;j>=1;--j){double b=cc[j]+2*x*b1-b2;b2=b1;b1=b;}
   value+=w[k]*(cc[0]+x*b1-b2);
  }out[i]=value;
 }
}

extern "C" void zstock_curve(int nt,int np,int nc,const double*x,const double*dw,
 const double*dt,const int*seg,const double*norm,const double*p,
 double rate,double sensitivity,uint64_t seed,double*out,double*terminal_z){
 const int nd=nc+8;
 std::vector<double> decay(nt),noise(nt),inverse_m(nt);
 for(int t=0;t<nt;++t){decay[t]=std::exp(-rate*dt[t]);noise[t]=std::sqrt(-std::expm1(-2*rate*dt[t]));inverse_m[t]=std::exp(-norm[8*t]);}
 auto worker=[&](int begin,int end){std::vector<double> logs(end-begin,0.),z(end-begin,0.);
  for(int t=0;t<nt;++t)for(int i=begin;i<end;++i){int j=i-begin;const double*state=x+3*(uint64_t(t)*np+i);
   double y=state[1]-p[2],f=state[2]-p[4],eta=p[0]*state[0]+p[1]*y*y+p[3]*std::exp(sensitivity*std::tanh(z[j]))*f*f,g,dg;
   hg(eta,p,g,dg);double variance=p[nd-1]+(p[7+seg[t]]-p[nd-1])*g*inverse_m[t];
   logs[j]+=-.5*variance*dt[t]+std::sqrt(variance)*dw[uint64_t(t)*np+i];
   if(sensitivity!=0)z[j]=decay[t]*z[j]+noise[t]*innovation(seed,uint64_t(t)*np+i);
  }
  for(int i=begin;i<end;++i){out[i*(nd+1)]=std::exp(logs[i-begin]);terminal_z[i]=z[i-begin];}
 };
 std::vector<std::thread> workers;for(int i=0;i<4;++i)workers.emplace_back(worker,np*i/4,np*(i+1)/4);for(auto&t:workers)t.join();
}
