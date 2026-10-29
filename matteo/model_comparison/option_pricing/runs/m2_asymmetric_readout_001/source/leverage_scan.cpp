// Fixed exogenous-factor leverage scan. Four sufficient statistics per grid.
#include "two_speed_precision_kernel.cpp"
static double sp(double x){return x>0?x+std::log1p(std::exp(-x)):std::log1p(std::exp(x));}
static double residual(double a,double h){double x=a*h;return h*(x<.01?x*x*(1./12-x*x/120+17*x*x*x*x/20160):1-std::tanh(x/2)/(x/2));}
extern "C" void leverage_stats(int nt,int pairs,int nd,int ne,const int*tr,const double*dec,const double*root,const double*cross,const double*mean,const double*proj,const double*time,const double*dt,const int*seg,const int*coarse,const int*reference,const int*expiry,const double*norm,const double*p,const double*x,uint64_t seed,double*out){
 Model m(x);std::vector<double>ds(nt),df(nt),ns(nt),nf(nt),cs(nt),cf(nt),rs(nt),rf(nt),qs(nt),qf(nt);
 for(int t=0;t<nt;++t){double h=dt[t];ds[t]=std::exp(-m.as*h);df[t]=std::exp(-m.af*h);ns[t]=std::sqrt(-std::expm1(-2*m.as*h));nf[t]=std::sqrt(-std::expm1(-2*m.af*h));rs[t]=std::sqrt(residual(m.as,h));rf[t]=std::sqrt(residual(m.af,h));cs[t]=std::sqrt(h-rs[t]*rs[t]);cf[t]=std::sqrt(h-rf[t]*rf[t]);qs[t]=-std::expm1(-2*m.as*time[t]);qf[t]=-std::expm1(-2*m.af*time[t]);}
 auto worker=[&](int begin,int end){for(int i=begin;i<end;++i){std::vector<double>s(nd,0),next(nd),z(nd+1);double zs=0,zf=0,sum[2][3][4]={},hold[2][3]={};
  for(int t=0;t<nt;++t){uint64_t k=(uint64_t(t)*pairs+i)*(nd+3);for(int j=0;j<nd+1;++j)z[j]=gaussian(seed,k+j);double dw=std::sqrt(dt[t])*z[0],es=gaussian(seed,k+nd+1),ef=gaussian(seed,k+nd+2);double bs=cs[t]*es+rs[t]*gaussian(seed^UINT64_C(0x2ad82d187dccb289),2*(uint64_t(t)*pairs+i)),bf=cf[t]*ef+rf[t]*gaussian(seed^UINT64_C(0x2ad82d187dccb289),2*(uint64_t(t)*pairs+i)+1);double pr[3]={};for(int j=0;j<nd;++j)for(int c=0;c<3;++c)pr[c]+=proj[c*nd+j]*s[j];
   for(int anti=0;anti<2;++anti){double sign=anti?-1:1,C=mean[3*t]+sign*pr[0],Y=mean[3*t+1]+sign*pr[1],F=mean[3*t+2]+sign*pr[2];double eta=p[0]*C+p[1]*std::pow(Y-p[2],2)+p[3]*std::pow(F-p[4],2);double old=std::exp(2*(p[5]-sp(p[6]*(p[5]-eta))/p[6])-norm[t]);double stress=(1-m.r)*std::exp(sign*m.bs*zs-.5*m.bs*m.bs*qs[t])+m.r*std::exp(sign*m.bf*zf-.5*m.bf*m.bf*qf[t]);double v=m.f+(.2-m.f)*x[seg[t]]*((1-m.w)*old+m.w*stress);
    if(coarse[t])hold[anti][0]=v;hold[anti][1]=v;if(reference[t])hold[anti][2]=m.f+(.2-m.f)*x[seg[t]]*((1-m.w)+m.w*stress);
    for(int g=0;g<3;++g){double vv=hold[anti][g],sv=std::sqrt(vv)*sign;sum[anti][g][0]+=vv*dt[t];sum[anti][g][1]+=sv*dw;sum[anti][g][2]+=sv*bs;sum[anti][g][3]+=sv*bf;}
    if(expiry[t]>=0){double*o=out+((uint64_t(anti*pairs+i)*ne+expiry[t])*12);for(int g=0;g<3;++g)for(int c=0;c<4;++c)o[g*4+c]=sum[anti][g][c];}
   }
   const double*r=root+uint64_t(tr[t])*nd*nd;for(int j=0;j<nd;++j){double v=cross[tr[t]*nd+j]*dw;for(int l=0;l<nd;++l)v+=r[j*nd+l]*z[l+1];next[j]=dec[tr[t]*nd+j]*s[j]+v;}s.swap(next);zs=ds[t]*zs+ns[t]*es;zf=df[t]*zf+nf[t]*ef;
  }
 }};std::vector<std::thread>threads;for(int j=0;j<4;++j)threads.emplace_back(worker,pairs*j/4,pairs*(j+1)/4);for(auto&t:threads)t.join();
}
// Independent cheap reference. W0 is integrated analytically at evaluation.
extern "C" void leverage_reference(int nt,int pairs,int ne,const double*time,const double*dt,const int*seg,const int*expiry,const double*x,uint64_t seed,double*out){
 Model m(x);auto worker=[&](int begin,int end){for(int i=begin;i<end;++i){double zs=0,zf=0,sum[2][4]={};for(int t=0;t<nt;++t){double h=dt[t],rs=std::sqrt(residual(m.as,h)),rf=std::sqrt(residual(m.af,h));uint64_t k=4*(uint64_t(t)*pairs+i);double es=gaussian(seed,k),ef=gaussian(seed,k+1),bs=std::sqrt(h-rs*rs)*es+rs*gaussian(seed,k+2),bf=std::sqrt(h-rf*rf)*ef+rf*gaussian(seed,k+3);for(int anti=0;anti<2;++anti){double sign=anti?-1:1,v=m.f+(.2-m.f)*x[seg[t]]*m.response(sign*zs,sign*zf,time[t]);sum[anti][0]+=v*h;sum[anti][2]+=std::sqrt(v)*sign*bs;sum[anti][3]+=std::sqrt(v)*sign*bf;if(expiry[t]>=0){double*o=out+(uint64_t(anti*pairs+i)*ne+expiry[t])*4;for(int j=0;j<4;++j)o[j]=sum[anti][j];}}zs=std::exp(-m.as*h)*zs+std::sqrt(-std::expm1(-2*m.as*h))*es;zf=std::exp(-m.af*h)*zf+std::sqrt(-std::expm1(-2*m.af*h))*ef;}}};std::vector<std::thread>threads;for(int j=0;j<4;++j)threads.emplace_back(worker,pairs*j/4,pairs*(j+1)/4);for(auto&t:threads)t.join();
}
