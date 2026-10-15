"""Reproduce Appendix C checks. Python 3 + numpy/scipy/sympy/matplotlib.
Run from any directory; results are written beside this script.
All square roots below denote the local chart connected to positive real s.
This is a research prototype, not a global pole-search implementation.
"""
from pathlib import Path
import os, json
os.environ.setdefault('MPLCONFIGDIR','/tmp/a1-matplotlib')
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
import sympy as sy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
m=.13957; M=m*m; a=4*M; mr=.77526; width=.1491; R=5.
sqrt=np.lib.scimath.sqrt

def nodes(n):
    x,w=leggauss(n); return (x+1)/2,w/2

def segment(lo,hi,n=180):
    t,w=nodes(n); theta=np.pi*t/2
    u=lo+(hi-lo)*np.sin(theta)**2
    du=(hi-lo)*np.pi*np.sin(theta)*np.cos(theta)
    return u,w*du

def endpoints(s,u):
    b=(sqrt(s)-m)**2; c=(sqrt(s)+m)**2
    B=sqrt(b-u)*sqrt(c-u)*sqrt(u-a)/sqrt(u)/2
    vp=(s+3*M-u)/2+B
    qp=sqrt(vp-a)
    h=2*m*((s-M)/2-u)/(sqrt(u)*qp)
    return qp,h

def toy_split(s,n=180,wrong=False,bend=0):
    b=(sqrt(s)-m)**2; u0=(s-M)/2; total=0j
    for lo,hi,sgn in [(a,u0,1),(u0,b,-1)]:
        # Optional analytic contour bend within each chart, endpoints fixed.
        t,w=nodes(n); th=np.pi*t/2; x=np.sin(th)**2; dx=np.pi*np.sin(th)*np.cos(th)
        u=lo+(hi-lo)*(x+1j*bend*x*(1-x))
        du=(hi-lo)*(1+1j*bend*(1-2*x))*dx
        qp,h=endpoints(s,u)
        qm=h if wrong else sgn*h
        total+=np.sum(w*du*2*(qp-qm))
    return total

def toy_reverse(s,n=240,mass=m):
    aa=4*mass**2; b=(sqrt(s)-mass)**2; c=(sqrt(s)+mass)**2
    u,w=segment(aa,b,n)
    return np.sum(w*sqrt(b-u)*sqrt(c-u)/sqrt(u))

def correction(s,n=180):
    b=(sqrt(s)-m)**2; u0=(s-M)/2
    u,w=segment(u0,b,n); qp,h=endpoints(s,u)
    return -4*np.sum(w*h)

# Threshold-factor-cancelled version of the repository's f1r, no normalization.
p0=np.sqrt(mr*mr/4-M)
rho0=2*p0/mr/(8*np.pi)
C=np.sqrt(1/R**2+p0**2)/p0

def reduced(v,q,sign=-1,eta=1.):
    ff=(v-a)/(4*p0*p0)*(1/R**2+p0*p0)/(1/R**2+(v-a)/4)
    mg=eta*mr*width*(q/sqrt(v)/(8*np.pi))/rho0*ff
    return C/sqrt(v-a+4/R**2)/(mr*mr-v+sign*1j*mg)

def W(s,u,v):
    A=sqrt(s); b=sqrt(u); c=sqrt(v)
    return (-4*M**3+4*M*s*s-4*M*M*A*b+4*M*A**3*b-4*M*M*A*c
        +4*M*A**3*c-9*M*M*b*c+8*M*A*A*b*c+A**4*b*c
        +14*M*A*b*b*c+2*A**3*b*b*c+9*M*b**3*c-A*A*b**3*c
        -4*A*b**4*c-2*b**5*c+14*M*A*b*c*c+2*A**3*b*c*c
        +12*M*b*b*c*c-6*A*b**3*c*c-4*b**4*c*c+9*M*b*c**3
        -A*A*b*c**3-6*A*b*b*c**3-5*b**3*c**3-4*A*b*c**4-4*b*b*c**4-2*b*c**5)

def kernel(s,u,v):
    return W(s,u,v)/(((sqrt(s)+sqrt(u))**2-M)*((sqrt(s)+sqrt(v))**2-M))

def interference(s,n=100,wrong=False,swap=False,eta=1.):
    b=(sqrt(s)-m)**2; u0=(s-M)/2; total=0j
    t,wt=nodes(n)
    for lo,hi,sgn in [(a,u0,1),(u0,b,-1)]:
        u,wu=segment(lo,hi,n); u=u[:,None]
        qp,h=endpoints(s,u); qm=h if wrong else sgn*h
        q=qm+(qp-qm)*t; v=a+q*q
        sign=1 if not swap else -1
        integr=2*q*(qp-qm)*reduced(u,sqrt(u-a),sign,eta)*reduced(v,q,-sign,eta)*kernel(s,u,v)
        total+=np.sum(wu[:,None]*wt*integr)
    return total/(2*np.pi*(8*np.pi)**2*s)

def interference_direct(s,n=180):
    b=(sqrt(s)-m)**2; u,wu=segment(a,b,n);u=u[:,None]
    qp,h=endpoints(s,u); vm=a+h*h; vp=a+qp*qp
    t,wt=nodes(n);v=vm+(vp-vm)*t
    integr=(vp-vm)*reduced(u,sqrt(u-a),1)*reduced(v,sqrt(v-a),-1)*kernel(s,u,v)
    return np.sum(wu[:,None]*wt*integr)/(2*np.pi*(8*np.pi)**2*s)

def run():
    u,s,MM=sy.symbols('u s M')
    lam=s*s+u*u+MM*MM-2*(s*u+s*MM+u*MM)
    identity=sy.factor((s-5*MM-u)**2/4-(u-4*MM)*lam/(4*u))
    assert sy.simplify(identity-MM*(s-MM-2*u)**2/u)==0
    result={'identity':str(identity),'toy':[],'full_unnormalized':[]}
    for ss in [.2,.5,1.5,1.5-.05j,1.5-.6j]:
        x=toy_reverse(ss); y=toy_split(ss); z=toy_split(ss,wrong=True); d=correction(ss)
        row=dict(s=str(ss),reverse=str(x),split=str(y),wrong=str(z),difference=str(z-x),correction=str(d),relative_agreement=float(abs(x-y)/abs(x)))
        result['toy'].append(row)
        assert abs(x-y)<2e-10*max(1,abs(x))
        assert abs(z-y-d)<1e-11
        assert abs(toy_split(ss,bend=.12)-y)<2e-10
    result['threshold']=[]
    for eps in [1e-3,1e-4,1e-5]:
        ss=9*M+eps; delta=(np.sqrt(ss)-m)**2-a
        result['threshold'].append([eps,float(toy_reverse(ss).real/delta**1.5)])
    result['massless']=float(toy_reverse(1.5,n=600,mass=0).real)
    for ss in [.5,1.5,1.5-.02j]:
        p=interference(ss,n=180); p2=interference(ss,n=260)
        other=interference(ss,n=260,swap=True); w=interference(ss,n=260,wrong=True)
        row=dict(s=str(ss),physical=str(p2),wrong=str(w),relative_order_error=float(abs(p2-other)/abs(p2)),relative_refinement=float(abs(p2-p)/abs(p2)),relative_mismatch=float(abs(w-p2)/abs(p2)))
        if isinstance(ss,float): row['direct']=str(interference_direct(ss,n=300))
        result['full_unnormalized'].append(row)
        assert abs(p2-other)<2e-8*max(abs(p2),1e-6)
        assert abs(p-p2)<2e-8*max(abs(p2),1e-6)
    ss=1.5-.02j; step=1e-5
    dx=(interference(ss+step,n=180)-interference(ss-step,n=180))/(2*step)
    dy=(interference(ss+1j*step,n=180)-interference(ss-1j*step,n=180))/(2j*step)
    result['CR_relative']=float(abs(dx-dy)/abs(dx))
    # Below rho-pi resonance, turning the width off is a nonsingular limit.
    result['width_off']=[dict(eta=e,mismatch=float(abs(interference(.3,eta=e,wrong=True)-interference(.3,eta=e)))) for e in [1,.1,.01,0]]
    result['global_chart_diagnostic']=[]
    for ss in [1.5-.1j,1.5-.3j,1.5-.6j,(1.209-.288j)**2]:
        low=interference(ss,n=180); high=interference(ss,n=260); swapped=interference(ss,n=260,swap=True)
        result['global_chart_diagnostic'].append(dict(s=str(ss),value=str(high),relative_refinement=float(abs(low-high)/abs(high)),relative_exchange=float(abs(high-swapped)/abs(high))))
    (ROOT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    sv=np.linspace(9*M+.001,2,140)
    phys=np.array([toy_reverse(x,n=90).real for x in sv]); wrong=np.array([toy_split(x,n=90,wrong=True).real for x in sv])
    fig,ax=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
    ax[0].plot(sv,phys,label='Physical / reversed order');ax[0].plot(sv,wrong,'--',label='One analytic lift in outer variable');ax[0].set(xlabel=r'$s$ [GeV$^2$]',ylabel=r'$X$ [GeV$^3$]');ax[0].legend(fontsize=8)
    ys=np.linspace(0,.7,100);z=1.5-1j*ys
    xp=np.array([toy_reverse(x,n=100) for x in z]);xw=np.array([toy_split(x,n=100,wrong=True) for x in z])
    ax[1].plot(ys,xp.real,label='Physical Re');ax[1].plot(ys,xp.imag,label='Physical Im');ax[1].plot(ys,xw.real,'--',label='Other Re');ax[1].plot(ys,xw.imag,'--',label='Other Im');ax[1].set(xlabel=r'$-\mathrm{Im}\,s$ at Re $s=1.5$ [GeV$^2$]',ylabel=r'$X$ [GeV$^3$]');ax[1].legend(fontsize=8,ncol=2)
    fig.savefig(ROOT/'figures/toy_continuation.pdf');fig.savefig(ROOT/'figures/toy_continuation.png',dpi=160);plt.close(fig)
    ss=1.5; b=(np.sqrt(ss)-m)**2;u0=(ss-M)/2;us=np.linspace(a+1e-7,b-1e-7,500);qp,h=endpoints(ss,us)
    fig,ax=plt.subplots(figsize=(6,3.3),layout='constrained');ax.plot(us,np.abs(h),label=r'Physical $q_-=|h|$');ax.plot(us,h,'--',label=r'Analytic lift $h$');ax.axvline(u0,color='gray',lw=.8);ax.set(xlabel=r'$u=\sigma_1$ [GeV$^2$]',ylabel=r'$q_-$ [GeV]');ax.legend();fig.savefig(ROOT/'figures/root_lift.pdf');plt.close(fig)
if __name__=='__main__':run()
