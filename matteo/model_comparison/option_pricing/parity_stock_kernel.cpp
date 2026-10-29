#include "controlled_readout_kernel.cpp"
extern "C" void stock_values(int nt,int np,const double*x,const double*dw,const double*dt,const int*seg,const double*norm,const double*p,double*out){
 auto worker=[&](int first,int last){
  for(int i=first;i<last;++i)out[i]=0;
  for(int t=0;t<nt;++t){double xi=p[7+seg[t]],floor=p[10],M=std::exp(norm[8*t]);
   for(int i=first;i<last;++i){double a[8];mult(x+3*(t*np+i),p,a);double v=floor+(xi-floor)*a[0]/M;out[i]+=-.5*v*dt[t]+std::sqrt(v)*dw[t*np+i];}
  }
  for(int i=first;i<last;++i)out[i]=std::exp(out[i]);
 };
 std::vector<std::thread> workers;for(int j=0;j<4;++j)workers.emplace_back(worker,np*j/4,np*(j+1)/4);for(auto&t:workers)t.join();
}
