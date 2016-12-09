# --------------------------------------------------------------------------- #
# ROMBOUTS et al. 2015 VEGETATION MODEL
# --------------------------------------------------------------------------- #

from PyDSTool import *
DSargs = args(name='Rombouts_Vegetation_model')

# --------------------------------------------------------------------------- #
# PARAMETERS
# --------------------------------------------------------------------------- #

pars = {'Ct':500,'Q0':342.5,'p':0.3,'av':0.1,'ag':0.4,'alphamax':0.85,'alphamin':0.25,
'Tal':263,'Tau':300,'B0':200,'B1':2.5,'Topt':283,'k':0.004,'gmma':0.01}
DSargs.pars = pars

# --------------------------------------------------------------------------- #
# INITIAL CONDITIONS
# --------------------------------------------------------------------------- #

icdict = {'T':280,'A':0.9}
DSargs.ics = icdict
DSargs.tdomain = [0,100]

# --------------------------------------------------------------------------- #
# EQUATIONS
# --------------------------------------------------------------------------- #

dTstr = '((1-alpha(T,A))*Q0-R0(T))/Ct'
dAstr = 'beta(T)*A*(1-A)-gmma*A'
DSargs.varspecs = {'T': dTstr, 'A': dAstr}
DSargs.fnspecs  = {'alpha': (['T','A'], '(1-p)*alpha0(T)+p*(av*A+ag*(1-A))'),
                   'beta': (['T'],'max(0,1-k*(T-Topt)**2)'),
                   'alpha0': (['T'],'max(alpha01(T),alpha02(T))'),
                   'alpha01': (['T'],'alphamax+min(0,(alphamin-alphamax)*(T-Tal)/(Tau-Tal))'),
                   'alpha02': (['T'],'alphamin'),
                   'R0':(['T'],'B0+B1*(T-Topt)')}
                               
# --------------------------------------------------------------------------- #
# SETUP GENERATOR
# --------------------------------------------------------------------------- #
ode = Generator.Vode_ODEsystem(DSargs)
ode.set(pars = {'gmma': 0.01} )
ode.set(ics = {'T': 290,'A':0.1} )

# --------------------------------------------------------------------------- #
# RUN FORWARD
# --------------------------------------------------------------------------- #

print('Computing curve...')
start = clock()
traj = ode.compute('test_traj')
print('done in %.3f seconds!' % (clock()-start))
pts = traj.sample()

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
PCargs.freepars = ['gmma']
PCargs.StepSize = 1e-4
PCargs.MaxNumPoints = 130
PCargs.MaxStepSize = 1e-1
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
PC.initpoint = 'EQ1:P1'
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
a=plt.figure()
PC.display(('gmma','A'),stability=True,figure=a,linewidth=2)
plt.ylabel(r'Vegetation Cover $A$',fontsize=15)
plt.xlabel(r'Vegetation death rate $\gamma$',fontsize=15)
plt.title('Rombouts Vegetation model')
plt.show()
a=plt.figure()
PC.display(('gmma','T'),stability=True,figure=a,linewidth=2)
plt.ylabel(r'Temperature $T$',fontsize=15)
plt.xlabel(r'Vegetation death rate $\gamma$',fontsize=15)
plt.title('Rombouts Vegetation model')
plt.show()