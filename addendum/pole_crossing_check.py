"""Locate the inner rho pole on the undeformed sector chart at s=1.5-0.3i."""
import numpy as np
from scipy.optimize import root
import continuation_checks as c

def den(z):
    q=z[0]+1j*z[1]
    inverse=1/c.reduced(c.a+q*q,q,-1)
    return [inverse.real,inverse.imag]
r=root(den,[.71,-.08]); assert r.success
q=r.x[0]+1j*r.x[1]
print('q pole:',q,'v pole:',c.a+q*q)
s=1.5-.3j;b=(c.sqrt(s)-c.m)**2;u0=(s-c.M)/2
for lo,hi,sgn in [(c.a,u0,1),(u0,b,-1)]:
    def eq(z):
        t,y=z;u=lo+(hi-lo)*t
        qp,h=c.endpoints(s,u);qm=sgn*h
        d=qm+y*(qp-qm)-q
        return [d.real,d.imag]
    r=root(eq,[.6,.7]); assert r.success and np.all(r.x>0) and np.all(r.x<1)
    print('sector',sgn,'linear outer fraction and inner fraction:',r.x,'residual:',np.linalg.norm(r.fun))
