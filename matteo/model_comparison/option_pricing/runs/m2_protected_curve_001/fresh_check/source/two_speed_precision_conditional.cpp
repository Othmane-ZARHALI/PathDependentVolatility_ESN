// Exact one-step Black-value martingale control for the discrete target stock.
// The included, frozen auxiliary control definitions are shared, not modified.
#include "two_speed_precision_kernel.cpp"
static inline double black_value(double logS,double K,double logK,double I){
 double sd=std::sqrt(I),d1=(logS-logK+.5*I)/sd;
 return std::exp(logS)*.5*std::erfc(-d1/std::sqrt(2.))-K*.5*std::erfc(-(d1-sd)/std::sqrt(2.));
}
extern "C" void precision_conditional(int nt,int np,int nq,const double*old,const double*dw,const double*eps,const double*time,const double*dt,const int*seg,const int*checkpoint,const double*x,const double*totalIV,const double*strikes,double*out){
 Model m(x);std::vector<double>ds(nt),df(nt),ns(nt),nf(nt),logK(nq);for(int t=0;t<nt;++t){ds[t]=std::exp(-m.as*dt[t]);df[t]=std::exp(-m.af*dt[t]);ns[t]=std::sqrt(-std::expm1(-2*m.as*dt[t]));nf[t]=std::sqrt(-std::expm1(-2*m.af*dt[t]));}for(int j=0;j<nq;++j)logK[j]=std::log(strikes[j]);
 auto worker=[&](int begin,int end){for(int i=begin;i<end;++i){double zs=0,zf=0,logS=0,prefixIV=0,qr=1;double*o=out+uint64_t(i)*(nq+1);
  for(int t=0;t<nt;++t){uint64_t k=uint64_t(t)*np+i;if(checkpoint[t])qr=m.response(zs,zf,time[t]);double q=m.response(zs,zf,time[t])+(1-m.w)*(old[k]-1),v=m.f+(.2-m.f)*x[seg[t]]*q,vr=m.f+(.2-m.f)*x[seg[t]]*qr,remaining=totalIV[i]-prefixIV;
   // Algebraically Ia = remaining_ref_after_this_step + actual_v*dt.
   // Form its difference directly, avoiding cumulative subtraction cancellation.
   double nextActual=remaining+(v-vr)*dt[t];
   for(int j=0;j<nq;++j)o[j]+=black_value(logS,strikes[j],logK[j],nextActual)-black_value(logS,strikes[j],logK[j],remaining);
   prefixIV+=vr*dt[t];logS+=-.5*v*dt[t]+std::sqrt(v)*dw[k];zs=ds[t]*zs+ns[t]*eps[2*k];zf=df[t]*zf+nf[t]*eps[2*k+1];
  }
  o[nq]=std::exp(logS);
 }};
 std::vector<std::thread>threads;for(int j=0;j<4;++j)threads.emplace_back(worker,np*j/4,np*(j+1)/4);for(auto&t:threads)t.join();
}
extern "C" void precision_bridge_integrals(int nt,int np,const double*normals,const double*time,const double*dt,const int*seg,const int*nodes,const int*left,const int*right,const double*coeff,const double*x,double*out){
 Model m(x);int dim=2*(nt-1);
 auto worker=[&](int begin,int end){for(int i=begin;i<end;++i){std::vector<double>zs(nt,0),zf(nt,0);const double*n=normals+uint64_t(i)*dim;
  for(int j=0;j<nt-1;++j){int k=nodes[j],l=left[j],r=right[j];const double*c=coeff+6*j;zs[k]=c[0]*zs[l]+c[1]*zs[r]+c[2]*n[2*j];zf[k]=c[3]*zf[l]+c[4]*zf[r]+c[5]*n[2*j+1];}
  double*o=out+uint64_t(i)*16;for(int t=0;t<nt;++t)o[seg[t]]+=m.response(zs[t],zf[t],time[t])*dt[t];
 }};
 std::vector<std::thread>threads;for(int j=0;j<4;++j)threads.emplace_back(worker,np*j/4,np*(j+1)/4);for(auto&t:threads)t.join();
}
