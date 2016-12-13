# --------------------------------------------------------------------------- #
# ROMBOUTS et al. 2015 VEGETATION MODEL
# --------------------------------------------------------------------------- #

from PyDSTool import *
DSargs = args(name='Stommel_Ocean_model')

# --------------------------------------------------------------------------- #
# PARAMETERS
# --------------------------------------------------------------------------- #

pars = { 'eta1': 3,
        'eta2': 0.5,
        'eta3': 0.3}
DSargs.pars = pars

# --------------------------------------------------------------------------- #
# INITIAL CONDITIONS
# --------------------------------------------------------------------------- #

icdict = {'T':1,'S':1}
DSargs.ics = icdict
DSargs.tdomain = [0,100]

# --------------------------------------------------------------------------- #
# EQUATIONS
# --------------------------------------------------------------------------- #

DSargs.varspecs = {'T': 'eta1-T*(1+abs(T-S))',
                   'S': 'eta2-S*(eta3+abs(T-S))'}
#DSargs.fnspecs = {'Jacobian': (['t','T','S'],
#                                """[[ -1-abs(T-S)-T*(T-S)/abs(T-S+0.01), (T**2-T*S)/abs(T-S+0.01) ],
#                                    [ -S*(T-S)/abs(T-S+0.01), -eta3-abs(T-S)+(T-S)/abs(T-S+0.01)]]""")}
                               
# --------------------------------------------------------------------------- #
# SETUP GENERATOR
# --------------------------------------------------------------------------- #
ode = Generator.Vode_ODEsystem(DSargs)
ode.set(pars = {'eta2': 0.01} )
ode.set(ics = {'T': 1,'S':1} )

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
ax2.plot(pts['t'], pts['S'], 'r', linewidth=2)
ax2.set_ylabel('Salinity',fontsize=15)
ax2.set_xlabel('Time',fontsize=15)
#%%
# --------------------------------------------------------------------------- #
# SETUP CONTINUATION CLASS
# --------------------------------------------------------------------------- #

PC = ContClass(ode)
PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['eta2']
PCargs.StepSize = 1e-4
PCargs.MaxNumPoints = 150
PCargs.MaxStepSize = 1e-1
PCargs.MinStepSize = 1e-8
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
PCargs.SaveEigen = True
PCargs.SaveJacobian = True
PC.newCurve(PCargs)

print('Computing curve...')
start = clock()
PC['EQ1'].forward()
print('done in %.3f seconds!' % (clock()-start))

#%%
# --------------------------------------------------------------------------- #
# PLOT BIFURCATION DIAGRAM
# --------------------------------------------------------------------------- #

PC['EQ1'].cleanLabels()
#PC['EQ2'].cleanLabels()
a=plt.figure()
PC.display(('eta2','S'),stability=True,figure=a,linewidth=2)
plt.ylabel(r'Salinity $S$',fontsize=15)
plt.xlabel(r'Parameter $\eta_2$',fontsize=15)
plt.title('Stommel ocean model')
plt.show()
a=plt.figure()
PC.display(('eta2','T'),stability=True,figure=a,linewidth=2)
plt.ylabel(r'Temperature $T$',fontsize=15)
plt.xlabel(r'Parameter $\eta_2$',fontsize=15)
plt.title('Stommel ocean model')
plt.show()