// Experimental independent normalized amplitude AFTER the frozen original cap.
#include "stochastic_vov_kernel.cpp"
extern "C" void poststock_curve(int nt,int np,int nc,const double*x,const double*dw,
 const double*dt,const int*seg,const double*norm,const double*p,
 double rate,double gamma,uint64_t seed,double*out,double*terminal_z){
 const int nd=nc+8;
 std::vector<double> decay(nt),noise(nt),clock(nt),inverse_m(nt);double time=0;
 for(int t=0;t<nt;++t){decay[t]=std::exp(-rate*dt[t]);noise[t]=std::sqrt(-std::expm1(-2*rate*dt[t]));clock[t]=-std::expm1(-2*rate*time);time+=dt[t];inverse_m[t]=std::exp(-norm[8*t]);}
 auto worker=[&](int begin,int end){std::vector<double> logs(end-begin,0.),z(end-begin,0.);
  for(int t=0;t<nt;++t)for(int i=begin;i<end;++i){int j=i-begin;const double*state=x+3*(uint64_t(t)*np+i);
   double eta=p[0]*state[0]+p[1]*std::pow(state[1]-p[2],2)+p[3]*std::pow(state[2]-p[4],2),g,dg;hg(eta,p,g,dg);
   double L=std::exp(gamma*z[j]-.5*gamma*gamma*clock[t]);
   double v=p[nd-1]+(p[7+seg[t]]-p[nd-1])*g*inverse_m[t]*L;
   logs[j]+=-.5*v*dt[t]+std::sqrt(v)*dw[uint64_t(t)*np+i];
   z[j]=decay[t]*z[j]+noise[t]*innovation(seed,uint64_t(t)*np+i);
  }
  for(int i=begin;i<end;++i){out[i*(nd+1)]=std::exp(logs[i-begin]);terminal_z[i]=z[i-begin];}
 };
 std::vector<std::thread> workers;for(int i=0;i<4;++i)workers.emplace_back(worker,np*i/4,np*(i+1)/4);for(auto&t:workers)t.join();
}
