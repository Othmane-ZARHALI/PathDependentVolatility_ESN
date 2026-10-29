// New readout only; all frozen Gaussian factor transitions are reused.
#include "leverage_scan.cpp"
#include <algorithm>
#include <stdexcept>
static double G(double eta,const double*p){return std::exp(2*(p[5]-sp(p[6]*(p[5]-eta))/p[6]));}
// Integrate independent C by a cubic Hermite convolution table, then Y/F GH.
// offset: [C nodes; Y/F pairs]. Caller audits separability of C from Y/F.
extern "C" void asymmetric_moments(int n,int nc,int nyf,int na,const double*means,const double*co,const double*cw,const double*yf,const double*yw,const double*p,double center,double scale,const double*alpha,double*out){
 double low=p[0]*means[0];for(int i=1;i<n;++i)low=std::min(low,p[0]*means[3*i]);double cmin=*std::min_element(co,co+nc),high=p[5]-p[0]*cmin+40/p[6],cap=std::exp(2*p[5]);
 const int N=4096;double step=(high-low)/N;std::vector<double>val(N+1),der(N+1);
 if(high>low)for(int j=0;j<=N;++j){double z=low+j*step;for(int c=0;c<nc;++c){double eta=z+p[0]*co[c],g=G(eta,p),sig=1/(1+std::exp(p[6]*(eta-p[5])));val[j]+=cw[c]*g;der[j]+=cw[c]*2*g*sig;}}
 auto worker=[&](int begin,int end){for(int i=begin;i<end;++i){double*o=out+uint64_t(i)*na*2;for(int a=0;a<na*2;++a)o[a]=0;for(int j=0;j<nyf;++j){double Y=means[3*i+1]+yf[2*j],F=means[3*i+2]+yf[2*j+1],z=p[0]*means[3*i]+p[1]*std::pow(Y-p[2],2)+p[3]*std::pow(F-p[4],2),v=cap;
   if(z<high&&step>0){int k=std::min(N-1,std::max(0,int((z-low)/step)));double t=(z-low-k*step)/step,delta=val[k+1]-val[k],c1=step*der[k],c2=3*delta-step*(2*der[k]+der[k+1]),c3=-2*delta+step*(der[k]+der[k+1]);v=val[k]+t*(c1+t*(c2+t*c3));}
   double l=sp(-(Y-center)/scale);for(int a=0;a<na;++a){double D=std::exp(alpha[a]*l);o[2*a]+=yw[j]*D;o[2*a+1]+=yw[j]*v*D;}
  }}};std::vector<std::thread>threads;for(int j=0;j<4;++j)threads.emplace_back(worker,n*j/4,n*(j+1)/4);for(auto&t:threads)t.join();
}
extern "C" void asymmetric_stock(int nt,int pairs,int nd,int ne,int na,const int*tr,const double*dec,const double*root,const double*cross,const double*mean,const double*proj,const double*time,const double*dt,const int*seg,const int*coarse,const int*reference,const int*expiry,const double*norm,const double*p,const double*models,const double*alpha,const double*normalizers,double scale,uint64_t seed,double*out){
 Model m(models);std::vector<double>ds(nt),df(nt),ns(nt),nf(nt),cs(nt),cf(nt),rs(nt),rf(nt),qs(nt),qf(nt);for(int t=0;t<nt;++t){double h=dt[t];ds[t]=std::exp(-m.as*h);df[t]=std::exp(-m.af*h);ns[t]=std::sqrt(-std::expm1(-2*m.as*h));nf[t]=std::sqrt(-std::expm1(-2*m.af*h));rs[t]=std::sqrt(residual(m.as,h));rf[t]=std::sqrt(residual(m.af,h));cs[t]=std::sqrt(h-rs[t]*rs[t]);cf[t]=std::sqrt(h-rf[t]*rf[t]);qs[t]=-std::expm1(-2*m.as*time[t]);qf[t]=-std::expm1(-2*m.af*time[t]);}
 auto worker=[&](int begin,int end){for(int i=begin;i<end;++i){std::vector<double>s(nd,0),next(nd),z(nd+1),logs(4*na,0),hold(2*na,0);double zs=0,zf=0,refv[2]={},reflog[2]={};
 for(int t=0;t<nt;++t){uint64_t k=(uint64_t(t)*pairs+i)*(nd+3);for(int j=0;j<nd+1;++j)z[j]=gaussian(seed,k+j);double dw=std::sqrt(dt[t])*z[0],es=gaussian(seed,k+nd+1),ef=gaussian(seed,k+nd+2),bs=cs[t]*es+rs[t]*gaussian(seed^UINT64_C(0x2ad82d187dccb289),2*(uint64_t(t)*pairs+i)),bf=cf[t]*ef+rf[t]*gaussian(seed^UINT64_C(0x2ad82d187dccb289),2*(uint64_t(t)*pairs+i)+1),shock=.6*dw-.8*(std::sqrt(1-m.r)*bs+std::sqrt(m.r)*bf),pr[3]={};for(int j=0;j<nd;++j)for(int c=0;c<3;++c)pr[c]+=proj[c*nd+j]*s[j];
 for(int anti=0;anti<2;++anti){double sign=anti?-1:1,C=mean[3*t]+sign*pr[0],Y=mean[3*t+1]+sign*pr[1],F=mean[3*t+2]+sign*pr[2],old=G(p[0]*C+p[1]*std::pow(Y-p[2],2)+p[3]*std::pow(F-p[4],2),p)*std::exp(-norm[t]);double stress=(1-m.r)*std::exp(sign*m.bs*zs-.5*m.bs*m.bs*qs[t])+m.r*std::exp(sign*m.bf*zf-.5*m.bf*m.bf*qf[t]),H=(1-m.w)*old+m.w*stress,l=sp(-sign*pr[1]/scale);
 for(int a=0;a<na;++a){const double*x=models+23*a;double q=H*std::exp(alpha[a]*l)/normalizers[t*na+a],v=m.f+(.2-m.f)*x[seg[t]]*q;if(coarse[t])hold[anti*na+a]=v;for(int grid=0;grid<2;++grid){double vv=grid?v:hold[anti*na+a];logs[(anti*na+a)*2+grid]+=-.5*vv*dt[t]+std::sqrt(vv)*sign*shock;}}
 if(reference[t])refv[anti]=m.f+(.2-m.f)*models[seg[t]]*((1-m.w)+m.w*stress);reflog[anti]+=-.5*refv[anti]*dt[t]+std::sqrt(refv[anti])*sign*shock;
 if(expiry[t]>=0){double*o=out+((uint64_t(anti*pairs+i)*ne+expiry[t])*(na+1)*2);for(int a=0;a<na;++a)for(int grid=0;grid<2;++grid)o[2*a+grid]=std::exp(logs[(anti*na+a)*2+grid]);o[2*na]=o[2*na+1]=std::exp(reflog[anti]);}
 }
 const double*r=root+uint64_t(tr[t])*nd*nd;for(int j=0;j<nd;++j){double v=cross[tr[t]*nd+j]*dw;for(int l=0;l<nd;++l)v+=r[j*nd+l]*z[l+1];next[j]=dec[tr[t]*nd+j]*s[j]+v;}s.swap(next);zs=ds[t]*zs+ns[t]*es;zf=df[t]*zf+nf[t]*ef;
 }}};std::vector<std::thread>threads;for(int j=0;j<4;++j)threads.emplace_back(worker,pairs*j/4,pairs*(j+1)/4);for(auto&t:threads)t.join();
}
