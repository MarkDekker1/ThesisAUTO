# --------------------------------------------------------------------------- #
# PREAMBLE
# --------------------------------------------------------------------------- #

import PyDSTool as dst
import numpy as np
from matplotlib import pyplot as plt

# --------------------------------------------------------------------------- #
# NAME
# --------------------------------------------------------------------------- #
DSargs = dst.args(name='Stommel model')

# --------------------------------------------------------------------------- #
# PARAMETERS
# --------------------------------------------------------------------------- #
DSargs.pars = { 'eta1': 3,
                'eta2': 1,
                'eta3': 0.3}
                
# --------------------------------------------------------------------------- #
# HELPER FUNCTIONS
# --------------------------------------------------------------------------- #

                                    
DSargs.fnspecs = {'Jacobian': (['t','T','S'],
                                """[[ -1-abs(T-S)-T*(T-S)/abs(T-S+0.01), (T**2-T*S)/abs(T-S+0.01) ],
                                    [ -S*(T-S)/abs(T-S+0.01), -eta3-abs(T-S)+(T-S)/abs(T-S+0.01)]]""")}

# --------------------------------------------------------------------------- #
# ODES
# --------------------------------------------------------------------------- #

DSargs.varspecs = {'T': 'eta1-T*(1+abs(T-S))',
                   'S': 'eta2-S*(eta3+abs(T-S))'}
                  
# --------------------------------------------------------------------------- #
# INITIAL CONDITIONS
# --------------------------------------------------------------------------- #
DSargs.ics      = {'T': 1.5, 'S': 1.5}

DSargs.xdomain = {'T': [0.3, 4], 'S': [0.3,4]}
DSargs.tdomain = [0,30]                         # set the range of integration.
ode  = dst.Generator.Vode_ODEsystem(DSargs)     # an instance of the 'Generator' class.


traj = ode.compute('test_traj')
pts = traj.sample()
#evs = traj.getEvents('event_x_a')

# figure 1 is the time evolution of the two variables
plt.figure(1)
plt.plot(pts['t'], pts['T'], 'b', linewidth=2)
plt.plot(pts['t'], pts['S'], 'r', linewidth=2)

# figure 2 is the phase plane
plt.figure(2,figsize=(10,6))
# phase plane tools are in the Toolbox module
from PyDSTool.Toolbox import phaseplane as pp

# plot vector field, using a scale exponent to ensure arrows are well spaced
# and sized
pp.plot_PP_vf(ode, 'S', 'T', scale_exp=-0.5)

# only one fixed point, hence [0] at end.
# n=4 uses three starting points in the domain to find any fixed points, to an
# accuracy of eps=1e-8.
fp_coord = pp.find_fixedpoints(ode, n=4, eps=1e-8)[0]
fp = pp.fixedpoint_2D(ode, dst.Point(fp_coord), eps=1e-8)

# n=3 uses three starting points in the domain to find nullcline parts, to an
# accuracy of eps=1e-8, and a maximum step for the solver of 0.1 units.
# The fixed point found is also provided to help locate the nullclines.
nulls_x, nulls_y = pp.find_nullclines(ode, 'S', 'T', n=3, eps=1e-8,subdomain={
                                    'S':DSargs.xdomain['S'],'T':DSargs.xdomain['T']},
                                      max_step=0.1, fps=[fp_coord])


# plot the fixed point
pp.plot_PP_fps(fp)

# plot the nullclines
plt.plot(nulls_x[:,0], nulls_x[:,1], 'b',linewidth=3)
plt.plot(nulls_y[:,0], nulls_y[:,1], 'g',linewidth=3)

# plot the trajectory
plt.plot(pts['S'], pts['T'], 'k', linewidth=1)

# plot the event points
#plt.plot(evs['x'], evs['y'], 'rs')

plt.axis('tight')
#plt.title('Phase plane')
plt.xlabel('Salinity',fontsize=15)
plt.ylabel('Temperature',fontsize=15)
plt.tick_params(axis='both',which='major',labelsize=15)

#%%
# --------------------------------------------------------------------------- #
# PLOT ONE CURVE
# --------------------------------------------------------------------------- #
traj = ode.compute('polarization')              # integrate ODE
pts  = traj.sample(dt=0.01)                     # Data for plotting
plt.plot(pts['t'], pts['T'])
plt.xlabel('time')                              # Axes labels
plt.ylabel('Temperature')                           # ...
#plt.ylim([0,65])                                # Range of the y axis
plt.title(ode.name)                             # Figure title from model name

# --------------------------------------------------------------------------- #
# PLOT MANY CURVES
# --------------------------------------------------------------------------- #
plt.figure()
plt.hold(True)                                  # Sequences of plot commands will not clear the existing figures
for eta2 in np.arange(0.5,1.5,0.1):
    for T0 in np.arange(0,3,1):
        ode.set(pars={'eta2':eta2})
        ode.set( ics = { 'T': T0 } )                # Initial condition
        # Trajectories are called pol0, pol1, ...
        # sample them on the fly to create Pointset tmp
        tmp = ode.compute('pol%3i' % eta2).sample()
        plt.plot(tmp['t'], tmp['T'])
plt.xlabel('time')
plt.ylabel('Temperature')
plt.title(ode.name + ' multi ICs')
plt.show()
#%%
# --------------------------------------------------------------------------- #
# SETUP CONTINUATION CLASS
# --------------------------------------------------------------------------- #
ode.set(pars = {'eta2': 0.1} )       # Lower bound of the control parameter 'i'
ode.set(ics =  {'T': 2,'S':2} )       # Close to one of the steady states present for i=-220

PC = dst.ContClass(ode)                 # Set up continuation class

PCargs = dst.args(name='EQ1', type='EP-C')  # 'EP-C' stands for Equilibrium Point Curve. The branch will be labeled 'EQ1'.
PCargs.freepars     = ['eta2']                      # control parameter(s) (it should be among those specified in DSargs.pars)
PCargs.MaxNumPoints = 150                        # The following 3 parameters are set after trial-and-error
PCargs.MaxStepSize  = 0.1
PCargs.MinStepSize  = 1e-7
PCargs.StepSize     = 2e-2
PCargs.LocBifPoints = 'LP'                       # detect limit points / saddle-node bifurcations
PCargs.SaveEigen    = True                       # to tell unstable from stable branches

# --------------------------------------------------------------------------- #
# CREATE BIFURCATION DIAGRAM
# --------------------------------------------------------------------------- #
plt.figure()
a=PC.newCurve(PCargs)
PC['EQ1'].forward()
data=PC['EQ1'].display(['eta2','T'], stability=True, figure=3,linewidth=2)        # stable and unstable branches as solid and dashed curves, resp.
plt.ylabel('Temperature',fontsize=15)
#plt.ylim([1,3.5])
plt.xlabel(r'$\eta_2$',fontsize=15)
#plt.xlim([0,2])
plt.title('')
plt.show()
#%%
# --------------------------------------------------------------------------- #
# DIAGRAM WITH POSITIONS OF LOCAL BIFURCATION POINTS
# --------------------------------------------------------------------------- #
PCargs = dst.args(name='SN1', type='LP-C')
PCargs.initpoint    = 'EQ1:LP1'
PCargs.freepars     = ['eta2', 'eta1']
PCargs.MaxStepSize  = 0.1
PCargs.MinStepSize  = 1e-8
PCargs.LocBifPoints = ['CP']
PCargs.MaxNumPoints = 150
PC.newCurve(PCargs)
PC['SN1'].forward()
PC['SN1'].backward()
PC['SN1'].display(['i','gca'], figure=4)