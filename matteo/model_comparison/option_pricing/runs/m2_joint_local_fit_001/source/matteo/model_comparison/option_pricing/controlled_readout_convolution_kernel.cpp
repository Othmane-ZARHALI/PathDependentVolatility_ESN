#include <cmath>
#include <algorithm>
#include <thread>
#include <vector>
#include <array>
#include <limits>
#include <atomic>
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
static void direct_expectation(int n,int q,const double*means,const double*offset,const double*w,const double*p,double*out){
 auto worker=[&](int first,int last){
  for(int i=first;i<last;++i){double*o=out+8*i;for(int j=0;j<8;++j)o[j]=0;
   for(int k=0;k<q;++k){double x[3],a[8];for(int j=0;j<3;++j)x[j]=means[3*i+j]+offset[3*k+j];mult(x,p,a);for(int j=0;j<8;++j)o[j]+=w[k]*a[j];}}
 };
 if(n<32){worker(0,n);return;}
 std::vector<std::thread> workers;for(int j=0;j<4;++j)workers.emplace_back(worker,n*j/4,n*(j+1)/4);for(auto& t:workers)t.join();
}

// Same tensor-GH32 rule, with the independent LINEAR C dimension convolved
// on a checked cubic-Hermite table. No change to the finite quadrature rule.
static std::atomic<long long> used(0),fallback(0),largest_table(0);
extern "C" void convolution_counts(long long*out){out[0]=used.load();out[1]=fallback.load();out[2]=largest_table.load();}
struct Knot {double z;std::array<double,4> v,d;};
static Knot convolution(double z,const double*co,const double*cw,const double*p){
 Knot t;t.z=z;t.v.fill(0);t.d.fill(0);
 for(int j=0;j<32;++j){
  double eta=z+p[0]*co[j],gap=p[6]*(p[5]-eta),u=std::exp(-std::abs(gap)),lp=std::log1p(u),s=gap>=0?1/(1+u):u/(1+u),h=(gap>=0?eta:p[5])-lp/p[6],g=std::exp(2*h),ge=2*g*s;
  double gee=g*(4*s*s-2*p[6]*s*(1-s)),gl=2*g*(1-s),gle=2*g*s*(1-s)*(2+p[6]),hk=(lp+std::abs(gap)*u/(1+u))/(p[6]*p[6]),gk=2*g*hk,gke=4*g*s*hk+2*g*gap*s*(1-s)/p[6];
  double v[4]={g,co[j]*ge,gl,gk},d[4]={ge,co[j]*gee,gle,gke};
  for(int k=0;k<4;++k){t.v[k]+=cw[j]*v[k];t.d[k]+=cw[j]*d[k];}
 }return t;
}
static void hermite(const Knot&a,const Knot&b,double z,double*v,double&dz){
 double h=b.z-a.z,t=(z-a.z)/h;
 for(int k=0;k<4;++k){double delta=b.v[k]-a.v[k],c1=h*a.d[k],c2=3*delta-h*(2*a.d[k]+b.d[k]),c3=-2*delta+h*(a.d[k]+b.d[k]);v[k]=a.v[k]+t*(c1+t*(c2+t*c3));if(k==0)dz=(c1+t*(2*c2+3*t*c3))/h;}
}
extern "C" void expectation(int n,int q,const double*means,const double*offset,const double*w,const double*p,double*out){
 if(n<32||q!=32768){direct_expectation(n,q,means,offset,w,p,out);return;}
 double co[32],cw[32]={},yw[1024]={};
 for(int j=0;j<32;++j){co[j]=offset[3*j*1024];for(int k=0;k<1024;++k){cw[j]+=w[j*1024+k];yw[k]+=w[j*1024+k];
  if(std::abs(offset[3*(j*1024+k)]-co[j])>1e-11||std::abs(offset[3*(j*1024+k)+1]-offset[3*k+1])>1e-11||std::abs(offset[3*(j*1024+k)+2]-offset[3*k+2])>1e-11){direct_expectation(n,q,means,offset,w,p,out);return;}
 }}
 double low=std::numeric_limits<double>::infinity(),cmin=co[0];for(int i=0;i<n;++i)low=std::min(low,p[0]*means[3*i]);for(double c:co)cmin=std::min(cmin,c);
 double high=p[5]-p[0]*cmin+40/p[6],cap=std::exp(2*p[5]);
 if(low>=high){for(int i=0;i<n;++i){double*o=out+8*i;for(int k=0;k<8;++k)o[k]=0;o[0]=cap;o[6]=2*cap;}return;}
 std::vector<Knot> table;for(int j=0;j<=64;++j)table.push_back(convolution(low+(high-low)*j/64,co,cw,p));bool pass=false;
 for(int iteration=0;iteration<12;++iteration){std::vector<Knot> refined;refined.push_back(table[0]);pass=true;
  for(size_t j=0;j+1<table.size();++j){double mid=(table[j].z+table[j+1].z)/2,pred[4],dz;hermite(table[j],table[j+1],mid,pred,dz);Knot exact=convolution(mid,co,cw,p);bool bad=false;
   for(int k=0;k<4;++k)if(std::abs(pred[k]-exact.v[k])>(k==0?1e-9:1e-8)*(std::abs(exact.v[k])+(k==0?1e-5:1.)))bad=true;
   if(std::abs(dz-exact.d[0])>1e-7*(std::abs(exact.d[0])+1.))bad=true;
   // Extra non-midpoint check catches Hermite cancellation at a midpoint.
   double zz=.75*table[j].z+.25*table[j+1].z;hermite(table[j],table[j+1],zz,pred,dz);Knot quarter=convolution(zz,co,cw,p);
   if(std::abs(pred[0]-quarter.v[0])>1e-9*(std::abs(quarter.v[0])+1e-5))bad=true;
   if(bad){pass=false;refined.push_back(exact);}refined.push_back(table[j+1]);
  }
  table.swap(refined);if(pass)break;if(table.size()>8193)break;
 }
 largest_table.store(std::max(largest_table.load(),(long long)table.size()));
 if(!pass){fallback.fetch_add(1);direct_expectation(n,q,means,offset,w,p,out);return;}used.fetch_add(1);
 auto worker=[&](int first,int last){
  for(int i=first;i<last;++i){double*o=out+8*i;for(int j=0;j<8;++j)o[j]=0;
   for(int k=0;k<1024;++k){double y=means[3*i+1]+offset[3*k+1]-p[2],f=means[3*i+2]+offset[3*k+2]-p[4],z=p[0]*means[3*i]+p[1]*y*y+p[3]*f*f,v[4],dz;
    if(z>=high){v[0]=cap;v[1]=0;v[2]=2*cap;v[3]=0;dz=0;}
    else{auto it=std::upper_bound(table.begin(),table.end(),z,[](double x,const Knot&t){return x<t.z;});size_t j=std::min<size_t>(std::max<size_t>(it-table.begin(),1),table.size()-1)-1;hermite(table[j],table[j+1],z,v,dz);}
    double a[8]={v[0],means[3*i]*dz+v[1],y*y*dz,-2*p[1]*y*dz,f*f*dz,-2*p[3]*f*dz,v[2],v[3]};for(int j=0;j<8;++j)o[j]+=yw[k]*a[j];
   }
  }
 };
 std::vector<std::thread> workers;for(int j=0;j<4;++j)workers.emplace_back(worker,n*j/4,n*(j+1)/4);for(auto&t:workers)t.join();
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
