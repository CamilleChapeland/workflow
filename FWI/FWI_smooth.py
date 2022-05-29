import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
import devito
from devito import configuration
#configuration['log-level'] = 'WARNING'
import argparse

#import warnings
#warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('starting_model', type=str)
parser.add_argument('rdir', type=str)
parser.add_argument('min_freq', type=float)
parser.add_argument('max_freq', type=float)
parser.add_argument('dx', type=float)
parser.add_argument('dz', type=float)
parser.add_argument('nshots', type=int)
args = parser.parse_args()

###################################################################
######################### MAKE PARAMETERS #########################
nshots = args.nshots  # Number of shots to create gradient from
#nshots = 2  # Number of shots to create gradient from
nreceivers = nshots  # Number of receiver locations per shot
fwi_iterations = 10  # Number of outer FWI iterations
min_freq = args.min_freq
max_freq = args.max_freq
dx = args.dx
dz = args.dz

##############################################################
######################### MAKE MODEL #########################
from examples.seismic import demo_model, Model, plot_velocity, plot_perturbation
# Define true and initial model
shape = (151, 201)  # Number of grid point (nx, nz)
spacing = (dx, dz)  # Grid spacing in m. The domain size is now 1km by 1km
origin = (0., 0.)  # Need origin to define relative source and receiver locations

# Make velocities
v = np.empty(shape, dtype=np.float32)
JMI_v = np.empty(shape, dtype=np.float32)
smooth_v = np.empty(shape, dtype=np.float32)

tru_vel = np.fromfile('Results/%s/truvel.bin'%args.rdir, dtype=np.float32)
tru_vel = np.reshape(tru_vel , shape)
v = tru_vel[:shape[0], :shape[1]]/1000

start_modelsm = np.fromfile('Results/%s/%s.bin'%(args.rdir,args.starting_model), dtype=np.float32)
start_modelsm = np.reshape(start_modelsm, shape)
smooth_v = start_modelsm[:shape[0], :shape[1]]/1000

# Make SeismicModel devito objects 
model = Model(vp=v, origin=origin, shape=shape, spacing=spacing, nbl=75, space_order=2, bcs="mask")
model0 = Model(vp=smooth_v, origin=origin, shape=shape, spacing=spacing, nbl=75, space_order=2, bcs="mask")

###############################################################
######################### MAKE SOURCE #########################
# Define source parameters
t0 = 0.
tn = 200. 
f0 = 0.2
source_type = "Ricker"
dt = model.critical_dt/1000
a = 1
fc = args.min_freq/1000
fs = 1/dt

############################################################################
######################### MAKE AQUISITION GEOMETRY #########################
# Define acquisition geometry: source

source_locations = np.empty((nshots, 2), dtype=np.float32)
source_locations[:, 1] = 0.
source_locations[:, 0] = np.linspace(0., model.domain_size[0], num=nshots)

# Define acquisition geometry: receivers
# Initialize receivers for synthetic and imaging data
rec_coordinates = np.empty((nreceivers, 2))
rec_coordinates[:, 0] = np.linspace(0, model.domain_size[0], num=nreceivers)
rec_coordinates[:, 1] = 0.

from buttersrc import AcquisitionGeometry
# Make geometry
geometry = AcquisitionGeometry(model, rec_coordinates, source_locations, 
                t0, tn, f0=f0, src_type="Butter", fs=fs, fc=fc)

##################################################################
######################### MAKE TRUE DATA #########################
# Compute synthetic data with forward operator 
from examples.seismic.acoustic import AcousticWaveSolver
solver = AcousticWaveSolver(model, geometry, space_order=4)

##                                                                       ##
##                                                                       ##
###########################################################################
######################### FWI WITH JMI SM #################################
###########################################################################
##                                                                       ##  
###########################################################################
######################### FWI GRADIENT PROCEDURES #########################
from devito import Eq, Operator

# Computes the residual between observed and synthetic data into the residual
def compute_residual(residual, dobs, dsyn):
    if residual.grid.distributor.is_parallel:
        # If we run with MPI, we have to compute the residual via an operator
        # First make sure we can take the difference and that receivers are at the 
        # same position
        assert np.allclose(dobs.coordinates.data[:], dsyn.coordinates.data)
        assert np.allclose(residual.coordinates.data[:], dsyn.coordinates.data)
        # Create a difference operator
        diff_eq = Eq(residual, dsyn.subs({dsyn.dimensions[-1]: residual.dimensions[-1]}) -
                               dobs.subs({dobs.dimensions[-1]: residual.dimensions[-1]}))
        Operator(diff_eq)()
    else:
        # A simple data difference is enough in serial
        residual.data[:] = dsyn.data[:] - dobs.data[:]
    
    return residual

# Create FWI gradient kernel
from devito import Function, TimeFunction, norm
from examples.seismic import Receiver

def fwi_gradient(vp_in):  
    # Create symbols to hold the gradient
    grad = Function(name="grad", grid=model.grid)
    # Create placeholders for the data residual and data
    residual = Receiver(name='residual', grid=model.grid,
                        time_range=geometry.time_axis, 
                        coordinates=geometry.rec_positions)
    d_obs = Receiver(name='d_obs', grid=model.grid,
                     time_range=geometry.time_axis, 
                     coordinates=geometry.rec_positions)
    d_syn = Receiver(name='d_syn', grid=model.grid,
                     time_range=geometry.time_axis, 
                     coordinates=geometry.rec_positions)
    objective = 0.
    for i in range(nshots):
        print(i)        
        # Update source location
        geometry.src_positions[0, :] = source_locations[i, :]
        # Generate synthetic data from true model
        _, _, _ = solver.forward(vp=model.vp, rec=d_obs)
        # Compute smooth data and full forward wavefield u0
        _, u0, _ = solver.forward(vp=vp_in, save=True, rec=d_syn)
        
        # Compute gradient from data residual and update objective function 
        compute_residual(residual, d_obs, d_syn)
        objective += .5*norm(residual)**2
        solver.gradient(rec=residual, u=u0, vp=vp_in, grad=grad)
    
    return objective, grad

########################################################################
######################### FWI INITIAL GRADIENT #########################
# Compute gradient of initial model
ff, update = fwi_gradient(model0.vp)
ff=np.rint(ff)

from devito import mmax
from examples.seismic import plot_image

# Show what the update does to the model
alpha = .5 / mmax(update)

###########################################################################
######################### FWI SOLUTION CONSTRAINT #########################
from devito import Min, Max
# Define bounding box constraints on the solution.
def update_with_box(vp, alpha, dm, vmin=1.5, vmax=3.5):
    """
    Apply gradient update in-place to vp with box constraint
    Notes:
    ------
    For more advanced algorithm, one will need to gather the non-distributed
    velocity array to apply constrains and such.
    """
    update = vp + alpha * dm 
    update_eq = Eq(vp, Max(Min(update, vmax), vmin))
    Operator(update_eq)()

########################################################################
######################### FWI GRADIENT DESCENT #########################
from devito import mmax
# Run FWI with gradient descent
history = np.zeros((fwi_iterations, 1))
for i in range(0, fwi_iterations):
    # Compute the functional value and gradient for the current
    # model estimate
    phi, direction = fwi_gradient(model0.vp)
    
    # Store the history of the functional values
    history[i] = phi
    
    # Artificial Step length for gradient descent
    # In practice this would be replaced by a Linesearch (Wolfe, ...)
    # that would guarantee functional decrease Phi(m-alpha g) <= epsilon Phi(m)
    # where epsilon is a minimum decrease constant
    alpha = .05 / mmax(direction)
    
    # Update the model estimate and enforce minimum/maximum values
    update_with_box(model0.vp , alpha , direction)
 
    # Log the progress made
    print('Objective value is %f at iteration %d' % (phi, i+1))

    np.save('Results/%s/%s-vel_model%d'%(args.rdir, args.starting_model, i+1), model0.vp.data)

################################################################
######################### PLOT RESULTS #########################
np.save('Results/%s/%s-final_vel_model'%(args.rdir, args.starting_model), model0.vp.data)

# Plot objective function decrease
plt.figure()
plt.loglog(history)
plt.xlabel('Iteration number')
plt.ylabel('Misift value Phi')
plt.title('Convergence')
plt.savefig('Results/%s/%slog_history.png'%(args.rdir, args.min_freq))
