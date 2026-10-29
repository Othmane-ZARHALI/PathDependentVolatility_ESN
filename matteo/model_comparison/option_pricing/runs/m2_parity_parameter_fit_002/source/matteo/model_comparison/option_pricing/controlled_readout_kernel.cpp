#include <cmath>
#include <algorithm>
#include <thread>
#include <vector>
// Columns: value, then derivatives wrt cC,qY,beta,qJ,betaJ,L,k.
static inline void mult(const double* x,const double* p,double* a){
 double y=x[1]-p[2],f=x[2]-p[4];
 double eta=p[0]*x[0]+p[1]*y*y+p[3]*f*f;
 double z=p[6]*(p[5]-eta);
 double u=std::exp(-std::abs(z)),lp=std::log1p(u);
 double s=z>=0?1./(1.+u):u/(1.+u);
 double h=(z>=0?eta:p[5])-lp/p[6];
 a[0]=std::exp(2*h);
 double d=2*a[0]*s;
 a[1]=d*x[0];a[2]=d*y*y;a[3]=-2*d*p[1]*y;
 a[4]=d*f*f;a[5]=-2*d*p[3]*f;
 a[6]=2*a[0]*(1-s);a[7]=2*a[0]*(lp+std::abs(z)*u/(1+u))/(p[6]*p[6]);
}
extern "C" void multiplier(int n,const double*x,const double*p,double*out){for(int i=0;i<n;++i)mult(x+3*i,p,out+8*i);}
extern "C" void expectation(int n,int q,const double*means,const double*offset,const double*w,const double*p,double*out){
 auto worker=[&](int first,int last){
  for(int i=first;i<last;++i){double*o=out+8*i;for(int j=0;j<8;++j)o[j]=0;
   for(int k=0;k<q;++k){double x[3],a[8];for(int j=0;j<3;++j)x[j]=means[3*i+j]+offset[3*k+j];mult(x,p,a);for(int j=0;j<8;++j)o[j]+=w[k]*a[j];}}
 };
 if(n<32){worker(0,n);return;}
 std::vector<std::thread> workers;for(int j=0;j<4;++j)workers.emplace_back(worker,n*j/4,n*(j+1)/4);for(auto& t:workers)t.join();
}
// norm columns: log M, dlogM for seven readout parameters.
extern "C" void stock(int nt,int np,int ne,const double*x,const double*dw,const double*dt,const int*seg,const int*events,const double*norm,const double*p,double*out){
 auto worker=[&](int first,int last){
 double*logs=new double[np]();double*grad=new double[np*11]();int ev=0;
 for(int t=0;t<nt;++t){int s=seg[t];double xi=p[7+s],floor=p[10],M=std::exp(norm[8*t]);
  for(int i=first;i<last;++i){double a[8];mult(x+3*(t*np+i),p,a);double ratio=a[0]/M,v=floor+(xi-floor)*ratio,sv=std::sqrt(v),shock=dw[t*np+i];
   logs[i]+=-.5*v*dt[t]+sv*shock;double fac=-.5*dt[t]+.5*shock/sv;
   for(int j=0;j<7;++j)grad[i*11+j]+=fac*(xi-floor)*(a[1+j]/M-ratio*norm[8*t+1+j]);
   grad[i*11+7+s]+=fac*ratio;grad[i*11+10]+=fac*(1-ratio);
  }
  if(ev<ne&&t+1==events[ev]){for(int i=first;i<last;++i){double*o=out+(ev*np+i)*12;o[0]=std::exp(logs[i]);for(int j=0;j<11;++j)o[1+j]=o[0]*grad[i*11+j];}++ev;}
 }
 delete[]logs;delete[]grad;
 };
 std::vector<std::thread> workers;for(int j=0;j<4;++j)workers.emplace_back(worker,np*j/4,np*(j+1)/4);for(auto& t:workers)t.join();
}
