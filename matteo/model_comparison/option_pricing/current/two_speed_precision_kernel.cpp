// Fixed-stress SPX control experiment. No fitted-model parameters are changed.
#include <cmath>
#include <vector>
#include <thread>
#include <cstdint>
static inline uint64_t hash64(uint64_t x){x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);return x^(x>>31);}
static inline double gaussian(uint64_t seed,uint64_t i){uint64_t a=hash64(seed+2*i*UINT64_C(0x9e3779b97f4a7c15)),b=hash64(seed+(2*i+1)*UINT64_C(0x9e3779b97f4a7c15));double u=(double(a>>11)+.5)*0x1.0p-53,v=(double(b>>11)+.5)*0x1.0p-53;return std::sqrt(-2*std::log(u))*std::cos(6.2831853071795864769*v);}
struct Model{
 double as,af,bs,bf,w,r,f;
 explicit Model(const double*x){double L=std::log(365.),lf=L*x[21],ls=lf+(L-lf)*x[18];as=std::log(2.)*365/std::exp(ls);af=std::log(2.)*365/std::exp(lf);bs=5*x[16];bf=5*x[20];w=x[17];r=x[22];f=std::pow(.05+.07*x[19],2);}
 double response(double zs,double zf,double t)const{return (1-w)+w*((1-r)*std::exp(bs*zs-.5*bs*bs*(-std::expm1(-2*as*t)))+r*std::exp(bf*zf-.5*bf*bf*(-std::expm1(-2*af*t))));}
};
extern "C" void precision_stock(int nt,int np,int nc,const double*old,const double*dw,const double*eps,const double*time,const double*dt,const int*seg,const int*checkpoint,const double*x,double*out,double*path){
 Model m(x);std::vector<double>ds(nt),df(nt),ns(nt),nf(nt);for(int t=0;t<nt;++t){ds[t]=std::exp(-m.as*dt[t]);df[t]=std::exp(-m.af*dt[t]);ns[t]=std::sqrt(-std::expm1(-2*m.as*dt[t]));nf[t]=std::sqrt(-std::expm1(-2*m.af*dt[t]));}
 auto worker=[&](int begin,int end){
  for(int i=begin;i<end;++i){double zs=0,zf=0,logS=0,logR=0,iv=0,qr=1;double gs[16]={},gr[16]={},gi[16]={};int c=0;
   for(int t=0;t<nt;++t){uint64_t k=uint64_t(t)*np+i;int j=seg[t];
    if(checkpoint[t]){qr=m.response(zs,zf,time[t]);uint64_t p=(uint64_t(c)*np+i)*3;path[p]=std::exp(logS);path[p+1]=std::exp(logR);path[p+2]=iv;++c;}
    double q=m.response(zs,zf,time[t])+(1-m.w)*(old[k]-1),v=m.f+(.2-m.f)*x[j]*q,vr=m.f+(.2-m.f)*x[j]*qr;
    double vs=std::sqrt(v),rs=std::sqrt(vr);logS+=-.5*v*dt[t]+vs*dw[k];logR+=-.5*vr*dt[t]+rs*dw[k];iv+=vr*dt[t];
    gs[j]+=(.2-m.f)*q*(-.5*dt[t]+.5*dw[k]/vs);gr[j]+=(.2-m.f)*qr*(-.5*dt[t]+.5*dw[k]/rs);gi[j]+=(.2-m.f)*qr*dt[t];
    zs=ds[t]*zs+ns[t]*eps[2*k];zf=df[t]*zf+nf[t]*eps[2*k+1];
   }
   double S=std::exp(logS),R=std::exp(logR);double*o=out+uint64_t(i)*51;o[0]=S;o[1]=R;o[2]=iv;for(int j=0;j<16;++j){o[3+j]=S*gs[j];o[19+j]=R*gr[j];o[35+j]=gi[j];}
   uint64_t p=(uint64_t(nc-1)*np+i)*3;path[p]=S;path[p+1]=R;path[p+2]=iv;
  }
 };
 std::vector<std::thread>threads;for(int j=0;j<4;++j)threads.emplace_back(worker,np*j/4,np*(j+1)/4);for(auto&t:threads)t.join();
}
// Mean-one response integrals H_j on exactly the same reference grid.
// Antithetic partners form one uncertainty unit; seeds distinguish iid batches.
extern "C" void precision_integrals(int nt,int pairs,const double*time,const double*dt,const int*seg,const double*x,uint64_t seed,double*out){
 Model m(x);std::vector<double>ds(nt),df(nt),ns(nt),nf(nt);for(int t=0;t<nt;++t){ds[t]=std::exp(-m.as*dt[t]);df[t]=std::exp(-m.af*dt[t]);ns[t]=std::sqrt(-std::expm1(-2*m.as*dt[t]));nf[t]=std::sqrt(-std::expm1(-2*m.af*dt[t]));}
 auto worker=[&](int begin,int end){for(int i=begin;i<end;++i){double zs=0,zf=0;for(int t=0;t<nt;++t){uint64_t k=uint64_t(t)*pairs+i;double*o=out+uint64_t(i)*32;o[seg[t]]+=m.response(zs,zf,time[t])*dt[t];o[16+seg[t]]+=m.response(-zs,-zf,time[t])*dt[t];zs=ds[t]*zs+ns[t]*gaussian(seed,2*k);zf=df[t]*zf+nf[t]*gaussian(seed,2*k+1);}}};
 std::vector<std::thread>threads;for(int j=0;j<4;++j)threads.emplace_back(worker,pairs*j/4,pairs*(j+1)/4);for(auto&t:threads)t.join();
}
