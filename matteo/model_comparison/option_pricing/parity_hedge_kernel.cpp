// Frozen candidate controls, one CALL convention across all strikes.
// Exact zero target mean: predictable piecewise-constant delta times increments
// of the adapted exponential-Euler stock. No future state in the delta.
#include "controlled_readout_kernel.cpp"
static inline double Phi(double x){return .5*std::erfc(-x/std::sqrt(2.));}
extern "C" void frozen_hedges(int nt,int np,int nk,const double*x,const double*dw,
 const double*dt,const int*seg,const int*rebalance,const double*remaining,
 const double*norm,const double*p,const double*strikes,double*out){
 auto worker=[&](int first,int last){
  std::vector<double> logs(np,0),previous(np,1),deltas(np*nk*2,0);
  for(int i=first;i<last;++i)for(int j=0;j<1+2*nk;++j)out[i*(1+2*nk)+j]=0;
  for(int t=0;t<nt;++t){
   double xi=p[7+seg[t]],floor=p[10],M=std::exp(norm[8*t]);
   for(int i=first;i<last;++i){
    double a[8];mult(x+3*(t*np+i),p,a);double v=floor+(xi-floor)*a[0]/M;
    if(rebalance[t]){
     double current=std::exp(logs[i]),change=current-previous[i];
     for(int j=0;j<nk;++j){
      double*o=out+i*(1+2*nk)+1+2*j;double*d=deltas.data()+(i*nk+j)*2;
      o[0]+=d[0]*change;o[1]+=d[1]*change;
      double w=remaining[t],w2=w*v/xi;
      d[0]=Phi((logs[i]-std::log(strikes[j])+.5*w)/std::sqrt(w));
      d[1]=Phi((logs[i]-std::log(strikes[j])+.5*w2)/std::sqrt(w2));
     }
     previous[i]=current;
    }
    logs[i]+=-.5*v*dt[t]+std::sqrt(v)*dw[t*np+i];
   }
  }
  for(int i=first;i<last;++i){
   double current=std::exp(logs[i]),change=current-previous[i];out[i*(1+2*nk)]=current;
   for(int j=0;j<nk;++j){double*o=out+i*(1+2*nk)+1+2*j;double*d=deltas.data()+(i*nk+j)*2;o[0]+=d[0]*change;o[1]+=d[1]*change;}
  }
 };
 std::vector<std::thread> workers;for(int j=0;j<4;++j)workers.emplace_back(worker,np*j/4,np*(j+1)/4);for(auto&t:workers)t.join();
}
