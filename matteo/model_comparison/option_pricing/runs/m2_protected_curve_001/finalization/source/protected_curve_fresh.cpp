// Fresh target-Q retained OU/equity/stress paths, paired models and stock grids.
#include "two_speed_precision_kernel.cpp"
static inline double softplus0(double z){return z>0?z+std::log1p(std::exp(-z)):std::log1p(std::exp(z));}
extern "C" void protected_fresh(int nt,int pairs,int nd,int ne,int ntr,const int*transition,const double*decay,const double*root,const double*cross,const double*mean,const double*projection,const double*time,const double*dt,const int*seg,const int*coarse,const int*reference,const int*expiry,const double*lognorm,const double*readout,const double*models,uint64_t seed,double*out){
 Model m(models);int np=2*pairs;std::vector<double>dS(nt),dF(nt),nS(nt),nF(nt),qS(nt),qF(nt);for(int t=0;t<nt;++t){dS[t]=std::exp(-m.as*dt[t]);dF[t]=std::exp(-m.af*dt[t]);nS[t]=std::sqrt(-std::expm1(-2*m.as*dt[t]));nF[t]=std::sqrt(-std::expm1(-2*m.af*dt[t]));qS[t]=-std::expm1(-2*m.as*time[t]);qF[t]=-std::expm1(-2*m.af*time[t]);}
 auto worker=[&](int begin,int end){for(int i=begin;i<end;++i){std::vector<double>state(nd,0),next(nd),normal(nd+1);double zs=0,zf=0,logs[2][2][2]={},logref[2][2]={},iv[2][2]={},vhold[2][2]={},qr[2]={1,1};
  for(int t=0;t<nt;++t){int tr=transition[t],j=seg[t];uint64_t counter=(uint64_t(t)*pairs+i)*(nd+3);for(int k=0;k<nd+1;++k)normal[k]=gaussian(seed,counter+k);double dw=std::sqrt(dt[t])*normal[0],es=gaussian(seed,counter+nd+1),ef=gaussian(seed,counter+nd+2);double projected[3]={};for(int k=0;k<nd;++k)for(int c=0;c<3;++c)projected[c]+=projection[c*nd+k]*state[k];
   for(int anti=0;anti<2;++anti){double sign=anti?-1:1,C=mean[3*t]+sign*projected[0],Y=mean[3*t+1]+sign*projected[1],F=mean[3*t+2]+sign*projected[2];double eta=readout[0]*C+readout[1]*std::pow(Y-readout[2],2)+readout[3]*std::pow(F-readout[4],2);double old=std::exp(2*(readout[5]-softplus0(readout[6]*(readout[5]-eta))/readout[6])-lognorm[t]);double A=std::exp(sign*m.bs*zs-.5*m.bs*m.bs*qS[t]),B=std::exp(sign*m.bf*zf-.5*m.bf*m.bf*qF[t]),stress=(1-m.r)*A+m.r*B;
    if(reference[t])qr[anti]=(1-m.w)+m.w*stress;
    double Q=(1-m.w)*old+m.w*stress;
    for(int model=0;model<2;++model){const double*x=models+23*model;double v=m.f+(.2-m.f)*x[j]*Q,vr=m.f+(.2-m.f)*x[j]*qr[anti];if(coarse[t])vhold[anti][model]=v;double vv[2]={vhold[anti][model],v};for(int level=0;level<2;++level)logs[anti][model][level]+=-.5*vv[level]*dt[t]+std::sqrt(vv[level])*sign*dw;logref[anti][model]+=-.5*vr*dt[t]+std::sqrt(vr)*sign*dw;iv[anti][model]+=vr*dt[t];
     if(expiry[t]>=0){uint64_t k=(((uint64_t(anti*pairs+i)*ne+expiry[t])*2+model)*4);out[k]=std::exp(logs[anti][model][0]);out[k+1]=std::exp(logs[anti][model][1]);out[k+2]=std::exp(logref[anti][model]);out[k+3]=iv[anti][model];}
    }
   }
   const double*r=root+uint64_t(tr)*nd*nd;for(int k=0;k<nd;++k){double innovation=cross[tr*nd+k]*dw;for(int l=0;l<nd;++l)innovation+=r[k*nd+l]*normal[l+1];next[k]=decay[tr*nd+k]*state[k]+innovation;}state.swap(next);zs=dS[t]*zs+nS[t]*es;zf=dF[t]*zf+nF[t]*ef;
  }
 }};
 std::vector<std::thread>threads;for(int j=0;j<4;++j)threads.emplace_back(worker,pairs*j/4,pairs*(j+1)/4);for(auto&t:threads)t.join();
}
