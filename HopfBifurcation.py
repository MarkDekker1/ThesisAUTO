# Exercise: Hopf bifurcation by the Van der Pol system
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

# Parameters and initial conditions
mu=-0.2
tmax=100.
dt=0.1
nt=np.int(tmax/dt)
x0=1
y0=0

# Functions
def dxdt(mu,x,y):
    return mu*(1.-y**2)*x-y

def dydt(x):
    return x

# Saving equilibria
maxequil=[]
minequil=[]
muvec=[]
xvecs=[]
yvecs=[]
for mu in np.arange(-2,2,0.01):
    

    # Numerical propagation (Semi-implicit Euler)
    xvec=[x0]
    yvec=[y0]
    tvec=[0.]
    for i in range(0,nt):
        x=xvec[i]+dxdt(mu,xvec[i],yvec[i])*dt
        y=yvec[i]+dydt(x)*dt
        
        tvec.append(i*dt)
        xvec.append(x)
        yvec.append(y)
    maxequil.append(max(xvec[800:]))
    minequil.append(min(xvec[800:]))
    muvec.append(mu)
    xvecs.append(xvec)
    yvecs.append(yvec)

# Plot
plt.figure(figsize=(8,4),dpi=1000)
plt.plot(muvec,maxequil,'-',linewidth=3)
plt.plot(muvec,minequil,'-',linewidth=3)
plt.ylabel(r'x$_{equil}$',fontsize=15)
plt.xlabel(r'$\mu$',fontsize=15)
plt.tick_params(axis='both',which='major',labelsize=15)

# Save plot
plt.savefig('HopfBifurcation.pdf',bbox_inches='tight')


# Plot
plt.figure(figsize=(8,4),dpi=1000)
for i in [0,100,200,300,399]:
    X=xvecs[i]
    Y=yvecs[i]
    plt.plot(X,Y,'-',linewidth=1)
plt.ylabel(r'Y',fontsize=15)
plt.xlabel(r'X',fontsize=15)
plt.tick_params(axis='both',which='major',labelsize=15)

# Save plot
plt.savefig('HopfBifurcation_phasespace.pdf',bbox_inches='tight')