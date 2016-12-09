# --------------------------------------------------------------------------- #
# ROMBOUTS et al. 2015 VEGETATION MODEL
# --------------------------------------------------------------------------- #

from PyDSTool import *
DSargs = args(name='Aleina_Vegetation_model')

# --------------------------------------------------------------------------- #
# PARAMETERS
# --------------------------------------------------------------------------- #

Le=2.501*10**6
cp=1000

pars = {'Le':1,'deltat':1./24.,'cp':1000,'cps':1000,'rhou':0.72,'rhol':1.200,
        'rhow':1000,'rhos':1800,'hu':10000,'hl':1000,'Zt':0.1,'Zd':4,'Z0':10,
        'n':0.4,'W0':5,'beta':2,'eu':0.25,'el':0.22,'ew':0.6,'es':0.85,'ae':0.25,
        'av':0.18,'S':0.5,'taut':180,'taud':360,'Fsun':340,'acl':0.6,'CD':0.008,
        'sw':0.18,'sa':0.46,'sfc':0.56,'Emax':0.0003,'Rmax':0.0002,'Ks':0.1,
        'betaw':12.7,'taua':3,'us':6,'f':0.5,'g0':0.8,'mu2':0.5,'mu1':0.1,
        'ds':0.005,'T0':283,'sigma':5.67*10**(-8)}
DSargs.pars = pars

# --------------------------------------------------------------------------- #
# INITIAL CONDITIONS
# --------------------------------------------------------------------------- #

icdict = {'Tu':1,'Tl':1,'Qu':1,'Ql':1,'W':1,'Tt':1,'Td':1,'St':1,'Sd':1,'b':1}
DSargs.ics = icdict
DSargs.tdomain = [0,2]

# --------------------------------------------------------------------------- #
# EQUATIONS
# --------------------------------------------------------------------------- #

dTustr = '(Le*DW(Qu,Tu)/deltat+rhol*cp*hl*DT(Tt,Tl)/deltat+(1-el)*(eu+ew*W/W0)*es*sigma*Tt**4+(eu+ew*W/W0)*el*sigma*Tl**4-2*(eu+ew*W/W0)*sigma*Tu**4)/(rhou*hu*cp)'
dTlstr = '(Qs(Tt,Tl)-rhol*cp*hl*DT(Tt,Tl)/deltat+el*es*sigma*Tt**4+(eu+ew*W/W0)*el*sigma*Tu**4-2*el*sigma*Tl**4)/(rhol*hl*cp)'
dTtstr = '((1-(b*av+(1-b)*ae))*Frad(W)-Qs(Tt,Tl)-es*sigma*Tt**4+el*sigma*Tl**4+(1-el)*es*(eu+ew*W/W0)*sigma*Tu**4-Le*Ei(St,Ql,Tl)*((1-b)+S*b)-Le*b*R(Sd,Ql)-rhos*cps*(Zt+Zd)*(Tt-Td)/taut)/(rhos*cps*Zt)'
dTdstr = '(-(Zt+Zd)*(Td-Tt)/taut-(Zd+Z0)*(Td-T0)/taud)/(Zd)'

dQustr = '(rhol*hl*Dq(Tt,Tl)/deltat-DW(Qu,Tu)/deltat)/(rhou*hu)'
dQlstr = '(((1-b)+S*b)*Ei(St,Ql,Tl)+b*R(Sd,Ql)-rhol*hl*Dq(Tt,Tl)/deltat)/(rhol*hl)'
dStstr = '(P(Tt,Tl,W)-Ei(St,Ql,Tl)*((1-b)+S*b)-L(St)-rhow*n*Zt*DI(St)/deltat)/(rhow*n*Zt)'
dSdstr = '(rhow*n*Zd*DI(St)/deltat+L(St)-b*R(Sd,Ql)-L(Sd)-rhow*n*Zd*DI(Sd)/deltat)/(rhow*n*Zd)'

dWstr = 'DW(Qu,Tu)/deltat-P(Tt,Tl,W)'
dbstr = 'g(St,Sd)*b*(1-b)-mu(Sd)*b'
t
DSargs.varspecs = {'Tu': dTustr,'Tl': dTlstr,'Tt': dTtstr,'Td': dTdstr,
                   'Qu': dQustr,'Ql': dQlstr,'St': dStstr,'Sd': dSdstr,
                   'W': dWstr,'b': dbstr}
DSargs.fnspecs  = {'DW': (['Qu','Tu'], 'rhou*hu*(Qu-qsat(Tu))'),
                    'qsat': (['T'],'6.1094*exp(17.625*T/(T+243.04))'),
                    'Ei': (['St','Ql','Tl'],'max(0,Emax*(St-sw)/(1-sw))*(qsat(Tl)-Ql)/qsat(Tl)'),
                    'R': (['Sd','Ql'],'min(max(0,Rmax*(Sd-sw)/(sa-sw)),Rmax)'),
                    'Frad': (['W'],'(1-alphaw(W))*Fsun'),
                    'alphaw': (['W'],'acl*W/W0'),
                    'L': (['Sv'],'Ks*(exp(betaw*(Sv-sfc))-1)/(exp(betaw*(1-sfc))-1)'),
                    'DI': (['Sv'],'max(0,Sv+ds-1)'),
                    'g': (['St','Sd'],'g0*gp(Sd)*gg(St)'),
                    'gp': (['Sd'],'min(max(0,(Sd-sw)/(sa-sw)),1)'),
                    'gg': (['St'],'1+tanh((St-sa)/ds)'),
                    'mu': (['Sd'],'(mu1-mu2)/2.*(1-tanh((Sd-sw)/ds))+mu1'),
                    'Dq': (['Tt','Tl'],'cp*DT(Tt,Tl)/(beta*Le)'),
                    'DT': (['Tt','Tl'],'10'),
                    'Qs': (['Tt','Tl'],'rhol*cp*CD*abs(us)*(Tt-Tl)'),
                    'P': (['Tt','Tl','W'],'min(1,(f*rhou*hl*Dq(Tt,Tl)/deltat+W/W0))*W/deltat'),
                    }
# --------------------------------------------------------------------------- #
# SETUP GENERATOR
# --------------------------------------------------------------------------- #
ode = Generator.Vode_ODEsystem(DSargs)
#%%
# --------------------------------------------------------------------------- #
# RUN FORWARD
# --------------------------------------------------------------------------- #

print('Computing curve...')
start = clock()
traj = ode.compute('test_traj')
print('done in %.3f seconds!' % (clock()-start))
pts = traj.sample()
#%%
fig, (ax1, ax2)= plt.subplots(2,1,figsize=(10,6),sharex=True)
ax1.plot(pts['t'], pts['T'], 'b', linewidth=2)
ax1.set_ylabel('Temperature',fontsize=15)
ax2.plot(pts['t'], pts['A'], 'r', linewidth=2)
ax2.set_ylabel('Vegetation cover',fontsize=15)
ax2.set_xlabel('Time',fontsize=15)
#%%
# --------------------------------------------------------------------------- #
# SETUP CONTINUATION CLASS
# --------------------------------------------------------------------------- #

PC = ContClass(ode)
PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['beta']
PCargs.StepSize = 1e-4
PCargs.MaxNumPoints = 250
PCargs.MaxStepSize = 1e-2
PCargs.MinStepSize = 1e-4
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
PCargs.SaveEigen = True
PCargs.SaveJacobian = True
PC.newCurve(PCargs)

print('Computing curve...')
start = clock()
PC['EQ1'].backward()
print('done in %.3f seconds!' % (clock()-start))

PCargs = args(name='EQ2', type='EP-C')
PC.initpoint = 'EQ1:P2'
PCargs.freepars = ['gmma']
PCargs.StepSize = 1e-4
PCargs.MaxNumPoints = 300
PCargs.MaxStepSize = 1e-1
PCargs.MinStepSize = 1e-4
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
PCargs.SaveEigen = True
PCargs.SaveJacobian = True
PC.newCurve(PCargs)

print('Computing curve...')
start = clock()
PC['EQ2'].forward()
print('done in %.3f seconds!' % (clock()-start))

#%%
# --------------------------------------------------------------------------- #
# PLOT BIFURCATION DIAGRAM
# --------------------------------------------------------------------------- #

PC['EQ1'].cleanLabels()
PC['EQ2'].cleanLabels()
a=plt.figure(figsize=(10,4))
PC.display(('gmma','A'),stability=True,figure=a)
plt.ylabel(r'Vegetation Cover $A$',fontsize=15)
plt.xlabel(r'Vegetation death rate $\gamma$',fontsize=15)
plt.title('Rombouts Vegetation model')
plt.show()
a=plt.figure(figsize=(10,4))
PC.display(('gmma','T'),stability=True,figure=a)
plt.ylabel(r'Vegetation Cover $A$',fontsize=15)
plt.xlabel(r'Vegetation death rate $\gamma$',fontsize=15)
plt.title('Rombouts Vegetation model')
plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Define the expression whose roots we want to find

a = 0.5
R = 1.6

def func(tau):
    return R - ((1.0 - np.exp(-tau))/(1.0 - np.exp(-a*tau))) 

# Plot it

tau = np.linspace(-0.5, 1.5, 201)

#plt.plot(tau, func(tau))
#plt.xlabel("tau")
#plt.ylabel("expression value")
#plt.grid()
#plt.show()

# Use the numerical solver to find the roots

tau_initial_guess = 0.5
tau_solution = fsolve(func, tau_initial_guess)

print "The solution is tau = %f" % tau_solution
print "at which the value of the expression is %f" % func(tau_solution)