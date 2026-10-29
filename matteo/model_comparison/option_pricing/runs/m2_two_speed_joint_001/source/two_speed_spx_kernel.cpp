// Adapted stock integration for the unchanged two-speed variance mixture.
#include <cmath>
#include <vector>
#include <thread>
#include <cstdint>
extern "C" void two_speed_stock(int nt,int np,int derivatives,const double*old,
 const double*dw,const double*eps,const double*time,const double*dt,const int*seg,
 const double*x,double*out){
 const int nd=23;const double L=std::log(365.),lhf=L*x[21],lhs=lhf+(L-lhf)*x[18];
 const double af=std::log(2.)*365/std::exp(lhf),as=std::log(2.)*365/std::exp(lhs);
 const double bs=5*x[16],bf=5*x[20],w=x[17],r=x[22],sigma=.05+.07*x[19],floor=sigma*sigma;
 std::vector<double> ds(nt),df(nt),ns(nt),nf(nt),qs(nt),qf(nt),hs(nt),hf(nt);
 for(int t=0;t<nt;++t){ds[t]=std::exp(-as*dt[t]);df[t]=std::exp(-af*dt[t]);ns[t]=std::sqrt(-std::expm1(-2*as*dt[t]));nf[t]=std::sqrt(-std::expm1(-2*af*dt[t]));qs[t]=-std::expm1(-2*as*time[t]);qf[t]=-std::expm1(-2*af*time[t]);hs[t]=-2*as*time[t]*std::exp(-2*as*time[t]);hf[t]=-2*af*time[t]*std::exp(-2*af*time[t]);}
 auto worker=[&](int begin,int end){
  int n=end-begin;std::vector<double> zS(n,0.),zF(n,0.),hS(n,0.),hF(n,0.),logs(n,0.),grad(n*nd,0.);
  for(int t=0;t<nt;++t)for(int i=begin;i<end;++i){int ii=i-begin;double*g=&grad[ii*nd];uint64_t k=uint64_t(t)*np+i;int j=seg[t];double A=std::exp(bs*zS[ii]-.5*bs*bs*qs[t]),B=std::exp(bf*zF[ii]-.5*bf*bf*qf[t]),stress=(1-r)*A+r*B,Q=(1-w)*old[k]+w*stress,c=(.2-floor)*x[j],v=floor+c*Q,sv=std::sqrt(v);
   logs[ii]+=-.5*v*dt[t]+sv*dw[k];
   if(derivatives){double factor=-.5*dt[t]+.5*dw[k]/sv,HHs=A*(bs*hS[ii]-.5*bs*bs*hs[t]),HHf=B*(bf*hF[ii]-.5*bf*bf*hf[t]);
    g[j]+=factor*(.2-floor)*Q;
    g[16]+=factor*c*5*w*(1-r)*A*(zS[ii]-bs*qs[t]);g[17]+=factor*c*(stress-old[k]);g[18]+=factor*c*w*(1-r)*HHs*(L-lhf);
    g[19]+=factor*2*sigma*.07*(1-x[j]*Q);g[20]+=factor*c*5*w*r*B*(zF[ii]-bf*qf[t]);g[21]+=factor*c*L*w*((1-r)*(1-x[18])*HHs+r*HHf);g[22]+=factor*c*w*(B-A);
   }
   double eS=eps[2*k],eF=eps[2*k+1];hS[ii]=ds[t]*hS[ii]+as*dt[t]*ds[t]*zS[ii]-as*dt[t]*ds[t]*ds[t]/ns[t]*eS;
   hF[ii]=df[t]*hF[ii]+af*dt[t]*df[t]*zF[ii]-af*dt[t]*df[t]*df[t]/nf[t]*eF;zS[ii]=ds[t]*zS[ii]+ns[t]*eS;zF[ii]=df[t]*zF[ii]+nf[t]*eF;
  }
  for(int i=begin;i<end;++i){int ii=i-begin;double S=std::exp(logs[ii]);out[i*(nd+1)]=S;for(int j=0;j<nd;++j)out[i*(nd+1)+j+1]=derivatives?S*grad[ii*nd+j]:0.;}
 };
 std::vector<std::thread> threads;for(int j=0;j<4;++j)threads.emplace_back(worker,np*j/4,np*(j+1)/4);for(auto&t:threads)t.join();
}
