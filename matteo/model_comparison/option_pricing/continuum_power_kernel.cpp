// G^c conditional expectations. configure(0,c); c in [0,1]. Value-only.
// Separate library: no changes to the archived M2 cap or pricing kernels.
#include <cmath>
#include <algorithm>
#include <vector>
#include <thread>
static int mode=0; static double strength=0;
extern "C" void configure(int m,double s){mode=m;strength=s;}
static inline double soft(double x){return std::max(x,0.)+std::log1p(std::exp(-std::abs(x)));}
static inline void hg(double eta,const double*p,double &g,double &dg){
 double gap=p[5]-eta,k=p[6],z=k*gap,s=z>=0?1/(1+std::exp(-z)):std::exp(z)/(1+std::exp(z));
 double h=std::min(eta,p[5])-std::log1p(std::exp(-std::abs(z)))/k,dh=s;
 if(mode==1){double d=2*std::log(2.)/k,r=std::hypot(gap,d);h=p[5]-.5*(gap+r);dh=.5*(1+gap/r);}
 if(mode==2){double u=soft(-z)/k;h+=strength*std::log1p(u);dh+=strength*(1-s)/(1+u);}
 g=std::exp(2*strength*h);dg=2*strength*g*dh;
}
static inline double value(const double*x,const double*p){
 double y=x[1]-p[2],f=x[2]-p[4],g,d;
 if(mode==3){double py[7];std::copy(p,p+7,py);py[5]=3.;hg(p[1]*y*y,py,g,d);double a=g;py[5]=2.5;hg(p[3]*f*f,py,g,d);return std::exp(2*p[0]*x[0])*a*g;}
 hg(p[0]*x[0]+p[1]*y*y+p[3]*f*f,p,g,d);return g;
}
extern "C" void multiplier(int n,const double*x,const double*p,double*out){for(int i=0;i<n;++i){out[8*i]=value(x+3*i,p);for(int j=1;j<8;++j)out[8*i+j]=0;}}
extern "C" void direct(int n,int q,const double*m,const double*o,const double*w,const double*p,double*out){
 for(int i=0;i<n;++i){double v=0;for(int j=0;j<q;++j){double x[3];for(int k=0;k<3;++k)x[k]=m[3*i+k]+o[3*j+k];v+=w[j]*value(x,p);}out[8*i]=v;for(int j=1;j<8;++j)out[8*i+j]=0;}
}
extern "C" void expectation(int n,int q,const double*m,const double*o,const double*w,const double*p,double*out){
 if(n<32||q!=32768){direct(n,q,m,o,w,p,out);return;}
 double co[32],cw[32]={},yw[1024]={};
 for(int j=0;j<32;++j){co[j]=o[3*j*1024];for(int k=0;k<1024;++k){cw[j]+=w[j*1024+k];yw[k]+=w[j*1024+k];
 if(std::abs(o[3*(j*1024+k)]-co[j])>1e-11||std::abs(o[3*(j*1024+k)+1]-o[3*k+1])>1e-11||std::abs(o[3*(j*1024+k)+2]-o[3*k+2])>1e-11){direct(n,q,m,o,w,p,out);return;}}}
 double low=1e100,high=-1e100,step=.01;
 if(mode!=3)for(int i=0;i<n;++i)for(int k=0;k<1024;++k){double y=m[3*i+1]+o[3*k+1]-p[2],f=m[3*i+2]+o[3*k+2]-p[4],z=p[0]*m[3*i]+p[1]*y*y+p[3]*f*f;low=std::min(low,z);high=std::max(high,z);}
 std::vector<double> v,d;
 if(mode!=3){int sz=int(std::ceil((high-low)/step))+2;v.resize(sz);d.resize(sz);for(int i=0;i<sz;++i)for(int j=0;j<32;++j){double g,dg;hg(low+i*step+p[0]*co[j],p,g,dg);v[i]+=cw[j]*g;d[i]+=cw[j]*dg;}}
 auto worker=[&](int a,int b){for(int i=a;i<b;++i){double result=0,cm=0;if(mode==3)for(int j=0;j<32;++j)cm+=cw[j]*std::exp(2*p[0]*(m[3*i]+co[j]));
 for(int k=0;k<1024;++k){double y=m[3*i+1]+o[3*k+1]-p[2],f=m[3*i+2]+o[3*k+2]-p[4],g;
 if(mode==3){double x[3]={0,m[3*i+1]+o[3*k+1],m[3*i+2]+o[3*k+2]};g=cm*value(x,p);}
 else{double z=p[0]*m[3*i]+p[1]*y*y+p[3]*f*f,t=(z-low)/step;int j=std::min(int(t),int(v.size())-2);t-=j;double delta=v[j+1]-v[j],c1=step*d[j],c2=3*delta-step*(2*d[j]+d[j+1]),c3=-2*delta+step*(d[j]+d[j+1]);g=v[j]+t*(c1+t*(c2+t*c3));}
 result+=yw[k]*g;}out[8*i]=result;for(int j=1;j<8;++j)out[8*i+j]=0;}};
 std::vector<std::thread> workers;for(int i=0;i<4;++i)workers.emplace_back(worker,n*i/4,n*(i+1)/4);for(auto&t:workers)t.join();
}
extern "C" void stock(int,int,int,const double*,const double*,const double*,const int*,const int*,const double*,const double*,double*){}
extern "C" void stock_curve(int nt,int np,int nc,int derivatives,const double*x,const double*dw,const double*dt,const int*seg,const double*norm,const double*p,double*out){
 const int nd=nc+8;auto worker=[&](int a,int b){std::vector<double> logs(np,0.);for(int t=0;t<nt;++t)for(int i=a;i<b;++i){double ratio=value(x+3*(t*np+i),p)/std::exp(norm[8*t]),v=p[nd-1]+(p[7+seg[t]]-p[nd-1])*ratio;logs[i]+=-.5*v*dt[t]+std::sqrt(v)*dw[t*np+i];}for(int i=a;i<b;++i)out[i*(nd+1)]=std::exp(logs[i]);};
 std::vector<std::thread> workers;for(int i=0;i<4;++i)workers.emplace_back(worker,np*i/4,np*(i+1)/4);for(auto&t:workers)t.join();
}
