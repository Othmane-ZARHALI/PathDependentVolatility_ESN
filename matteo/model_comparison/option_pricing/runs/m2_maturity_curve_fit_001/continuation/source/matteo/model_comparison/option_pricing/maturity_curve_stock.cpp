#include "controlled_readout_kernel.cpp"
// Arbitrary deterministic curve on existing left-point stock integration grid.
extern "C" void stock_curve(int nt,int np,int nc,int derivatives,const double*x,const double*dw,const double*dt,const int*seg,const double*norm,const double*p,double*out){
 const int nd=nc+8,fi=nd-1;
 auto worker=[&](int first,int last){
  std::vector<double> logs(np,0.),grad(derivatives?np*nd:0,0.);
  for(int t=0;t<nt;++t){int s=seg[t];double xi=p[7+s],floor=p[fi],M=std::exp(norm[8*t]);
   for(int i=first;i<last;++i){double a[8];mult(x+3*(t*np+i),p,a);double ratio=a[0]/M,v=floor+(xi-floor)*ratio,sv=std::sqrt(v),shock=dw[t*np+i];
    logs[i]+=-.5*v*dt[t]+sv*shock;
    if(derivatives){double fac=-.5*dt[t]+.5*shock/sv;
     for(int j=0;j<7;++j)grad[i*nd+j]+=fac*(xi-floor)*(a[1+j]/M-ratio*norm[8*t+1+j]);
     grad[i*nd+7+s]+=fac*ratio;grad[i*nd+fi]+=fac*(1-ratio);
    }
   }
  }
  for(int i=first;i<last;++i){double*o=out+i*(nd+1);o[0]=std::exp(logs[i]);if(derivatives)for(int j=0;j<nd;++j)o[1+j]=o[0]*grad[i*nd+j];}
 };
 std::vector<std::thread> workers;for(int j=0;j<4;++j)workers.emplace_back(worker,np*j/4,np*(j+1)/4);for(auto&t:workers)t.join();
}
