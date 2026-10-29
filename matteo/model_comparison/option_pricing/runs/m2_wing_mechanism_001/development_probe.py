import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import norm

def EY(s,w=.1):
 a,b=.3*s,3*s
 f=lambda z:np.sqrt((1-w)*np.exp(a*z-a*a/2)+w*np.exp(b*z-b*b/2))*norm.pdf(z)
 return quad(f,-12,12,epsabs=1e-12)[0]
s=0.6
m=EY(s);a,b=.3*s,3*s
print(s,a,b,m)
for k in [1,2,4,8]:
 root=brentq(lambda z:np.sqrt(.9*np.exp(a*z-a*a/2)+.1*np.exp(b*z-b*b/2))/m-k,-12,12)
 c=quad(lambda z:(np.sqrt(.9*np.exp(a*z-a*a/2)+.1*np.exp(b*z-b*b/2))/m-k)*norm.pdf(z),root,12,epsabs=1e-12)[0]
 def black(v):
  d=-np.log(k)/v+v/2
  return ndtr(d)-k*ndtr(d-v)
 iv=brentq(lambda v:black(v)-c,.0001,10)
 print(k,c,iv)
print('single',np.sqrt(np.log(1.5)))
