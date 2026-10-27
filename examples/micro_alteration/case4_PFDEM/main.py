# this code presents a Phase-Field Discrete Element Method applied for a single grain under the pressure-solution phenomenon.
# last update: 09/09/2026
# authors: Alexandre Sac-Morane, alexandre.sac-morane@enpc.fr

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Librairies
# ------------------------------------------------------------------------------------------------------------------------------------------ #

import numpy as np
import matplotlib.pyplot as plt
import math, shutil, os
from pathlib import Path

# own librairies 
from parameters import get_parameters
from prepare_pf_lib import *
from run_pf_lib import *
from pp_lib import *

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------------------------------------------------------------ #

def define_ic(dict_user, dict_sample):
    '''
    Define the initial conditions for the simulation.

    dict_sample['x_L', 'y_L', 'pos_1', 'eta_1_map', 'c_map'] are created
    '''
    # create the mesh
    dict_sample['x_L'] = np.linspace(dict_user['x_min'], dict_user['x_max'], dict_user['n_mesh_x'])
    dict_sample['y_L'] = np.linspace(dict_user['y_min'], dict_user['y_max'], dict_user['n_mesh_y'])

    # position of the grain
    dict_sample['pos_1'] = [0, dict_user['radius']]

    # initialize and compute the phase map
    dict_sample['eta_1_map'] = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))
    # iterate over x
    for i_x in range(len(dict_sample['x_L'])):
        x = dict_sample['x_L'][i_x]
        # iterate over y
        for i_y in range(len(dict_sample['y_L'])):
            y = dict_sample['y_L'][i_y]
            # distance to grain center
            d_node_to_g1 = np.linalg.norm(np.array([x,y])-np.array(dict_sample['pos_1']))
            # compute the phase variable
            if d_node_to_g1 <= dict_user['radius']-dict_user['w_int']/2 : # inside the grain
                dict_sample['eta_1_map'][-1-i_y, i_x] = 1
            elif dict_user['radius']-dict_user['w_int']/2 < d_node_to_g1 and d_node_to_g1 < dict_user['radius']+dict_user['w_int']/2: # in the interface
                dict_sample['eta_1_map'][-1-i_y, i_x] = 0.5*(1+math.cos(math.pi*(d_node_to_g1-dict_user['radius']+dict_user['w_int']/2)/dict_user['w_int'])) # a cosine profile is assumed
            elif dict_user['radius']+dict_user['w_int']/2 <= d_node_to_g1 : # outside the grain
                dict_sample['eta_1_map'][-1-i_y, i_x] = 0 

    # initialize and compute the solute map
    dict_sample['c_map'] = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))
    for i_x in range(len(dict_sample['x_L'])):
        for i_y in range(len(dict_sample['y_L'])):
            dict_sample['c_map'][-1-i_y, i_x] = 1 # system at the equilibrium initialy

# ---------------------------------------------------------------------

def prepare_dem(dict_user, dict_sample):
    '''
    Prepare the DEM simulation.

    The shape of the grains is interpolated from the Phase-Field map.
    This shape can be represented by a polyhedral or by a Level-Set function.

    In the configuration investigated here, this step is not necessary as the home made solver (see run_dem()) uses directly the phase-field map.
    '''
    pass

# ---------------------------------------------------------------------

def run_dem(dict_user, dict_sample):
    '''
    Run the DEM simulation.

    Here an home made solver is used. More complex configurations could require the use of more developped solvers such as YADE.
    The idea is to determine the y coordinates that respect the imposed contact volume.
    '''
    # start the investigation from the bottom
    i_y = 0
    # compute the contact volume at this current coordinate
    contact_volume = dict_sample['eta_1_map'][-1-i_y,:].sum()*(dict_sample['x_L'][1]-dict_sample['x_L'][0])*(dict_sample['y_L'][1]-dict_sample['y_L'][0])*1
    # save the previous contact volume
    contact_volume_prev = 0
    # compare the contact volume to the imposed solicitation
    while contact_volume < dict_user['imposed_contact_volume']:
        # iterate over the y coordinates
        i_y = i_y + 1
        # update the contact volume and the previous contact volume
        contact_volume_prev = contact_volume
        contact_volume = contact_volume + dict_sample['eta_1_map'][-1-i_y,:].sum()*(dict_sample['x_L'][1]-dict_sample['x_L'][0])*(dict_sample['y_L'][1]-dict_sample['y_L'][0])*1
    # interpolate the y coordinate to respect the imposed contact volume
    # this value is the displacement of the grain computed by the DEM simulation
    dict_sample['displacement'] = np.array([0, -(dict_sample['y_L'][i_y-1] + (dict_user['imposed_contact_volume']-contact_volume_prev)/(contact_volume-contact_volume_prev)*(dict_sample['y_L'][i_y]-dict_sample['y_L'][i_y-1]))])

# ---------------------------------------------------------------------

def prepare_pf(dict_user, dict_sample):
    '''
    Prepare the Phase-Field simulation.

    The rigid body motion of the grain is applied to the phase map.
    The solute map is updated to ensure that no solute is located in the solid.
    The geometry of the contact is characterized (a focus is made on the contact surface).
    The solid activity map is computed.
    The diffusivity map is computed.
    '''
    # Update the phase map of the grain
    update_phase_map(dict_user, dict_sample) # see prepare_pf_lib.py

    # Update the solute map
    update_solute_map(dict_user, dict_sample) # see prepare_pf_lib.py

    # Characterize the contact
    characterize_contact(dict_user, dict_sample) # see prepare_pf_lib.py

    # Compute the solid activity map
    compute_as_map(dict_user, dict_sample) # see prepare_pf_lib.py

    # Compute the diffusivity map
    compute_diffusivity_map(dict_user, dict_sample) # see prepare_pf_lib.py

# ---------------------------------------------------------------------

def run_pf(dict_user, dict_sample):
    '''
    Run the Phase-Field simulation with the solver MOOSE.
    '''
    # create a folder dedicated to the MOOSE simulation
    create_folder('MOOSE_simulation')

    # write the variables maps required for the MOOSE simulation
    write_map_txt(dict_user, dict_sample, 'eta_1_map') # see run_pf_lib.py
    write_map_txt(dict_user, dict_sample, 'c_map') # see run_pf_lib.py
    write_map_txt(dict_user, dict_sample, 'as_map') # see run_pf_lib.py
    write_map_txt(dict_user, dict_sample, 'kc_map') # see run_pf_lib.py

    # write the MOOSE input file
    write_moose_input_file(dict_user, dict_sample) # see run_pf_lib.py

    # run MOOSE to solve the PF
    os.system('mpiexec -n '+str(dict_user['n_proc'])+' ~/projects/moose/modules/phase_field/phase_field-opt -i Case4_PFDEM.i')

    # sort the output
    os.rename('Case4_PFDEM.i', 'MOOSE_simulation/Case4_PFDEM.i')
    os.rename('Case4_PFDEM_out.e', 'MOOSE_simulation/Case4_PFDEM_out.e')
    # iterate on the.vtu (one per processor)
    for i_proc in range(dict_user['n_proc']):
        # configuration before PF simulation
        os.rename('Case4_PFDEM_other_000_'+str(i_proc)+'.vtu',\
                  'vtk/Case4_PFDEM_other_'+index_to_3str((dict_sample['i_PFDEM_ite']-1)*2)+'_'+str(i_proc)+'.vtu')
        # configuration after PF simulation
        os.rename('Case4_PFDEM_other_001_'+str(i_proc)+'.vtu',\
                    'vtk/Case4_PFDEM_other_'+index_to_3str((dict_sample['i_PFDEM_ite']-1)*2+1)+'_'+str(i_proc)+'.vtu')
    # write a new pvtu from the previous one
    adapt_pvtu('Case4_PFDEM_other_000.pvtu', index_to_3str((dict_sample['i_PFDEM_ite']-1)*2)) # see run_pf_lib.py
    adapt_pvtu('Case4_PFDEM_other_001.pvtu', index_to_3str((dict_sample['i_PFDEM_ite']-1)*2+1)) # see run_pf_lib.py
    # move the previous pvtu into the bin folder
    os.rename('Case4_PFDEM_other_000.pvtu',\
              'MOOSE_simulation/Case4_PFDEM_other_000.pvtu') 
    os.rename('Case4_PFDEM_other_001.pvtu',\
              'MOOSE_simulation/Case4_PFDEM_other_001.pvtu') 
    # read the maps
    read_vtk(dict_user, dict_sample, index_to_3str((dict_sample['i_PFDEM_ite']-1)*2+1))
    # clean the unused output
    shutil.rmtree('MOOSE_simulation')

# ---------------------------------------------------------------------

def create_folder(folder_name):
    '''
    Create a folder (if it already exists, it is earased).
    '''
    if Path(folder_name).exists():
        shutil.rmtree(folder_name)
    os.mkdir(folder_name)

# ---------------------------------------------------------------------

def index_to_3str(index):
    '''
    Convert an index into a str with the form 'xxx'.
    '''
    if index < 10:
        return '00'+str(index)
    if index < 100:
        return '0'+str(index)
    else :
        return str(index)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Generation of the initial conditions
# ------------------------------------------------------------------------------------------------------------------------------------------ #

# load the parameters from the dedicated file
dict_user, dict_sample = get_parameters()

# define the initial conditions
define_ic(dict_user, dict_sample)

# create folder for the vtk files
create_folder('vtk')

# initial pp
compute_distribution_mass(dict_user, dict_sample) # see pp_lib.py

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------------------------------------------------------------ #

# iterate over the PFDEM iterations
while dict_sample['i_PFDEM_ite'] < dict_user['n_PFDEM_ite']:
    dict_sample['i_PFDEM_ite'] = dict_sample['i_PFDEM_ite'] + 1
    # user interface
    print('\nStep',dict_sample['i_PFDEM_ite'],'/',dict_user['n_PFDEM_ite'],'\n')

    # prepare the DEM simulation
    prepare_dem(dict_user, dict_sample)

    # run the DEM simulation
    run_dem(dict_user, dict_sample)

    # prepare the Phase-Field simulation
    prepare_pf(dict_user, dict_sample)

    # run the Phase-Field simulation
    run_pf(dict_user, dict_sample)

    # post-proccess
    compute_cog(dict_user, dict_sample) # see pp_lib.py
    compute_as_max(dict_user, dict_sample) # see pp_lib.py
    compute_distribution_mass(dict_user, dict_sample) # see pp_lib.py

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Print the results
# ------------------------------------------------------------------------------------------------------------------------------------------ #

# cog
fig, ax1 = plt.subplots(1,1, figsize=(16,9))
ax1.plot(dict_sample['L_cog_y'])
ax1.set_xlabel('iteration (-)')
ax1.set_ylabel('coordinate y of the cog (-)')
fig.tight_layout()
fig.savefig('cog.png')
plt.close()

# as
fig, ax1 = plt.subplots(1,1, figsize=(16,9))
ax1.plot(dict_sample['L_as_max'])
ax1.set_xlabel('iteration (-)')
ax1.set_ylabel('solid activity at the contact (-)')
fig.tight_layout()
fig.savefig('as.png')
plt.close()

# mass
fig, ax1 = plt.subplots(1,1, figsize=(16,9))
ax1.plot(dict_sample['L_s_eta_1'], color='r', label='eta')
ax1.set_ylabel('sum eta (-)', color='r')
ax1b = ax1.twinx()
ax1b.plot(dict_sample['L_s_c'], color='b', label='c')
ax1b.set_ylabel('sum c (-)', color='b')
ax1.set_xlabel('iteration (-)')
fig.tight_layout()
fig.savefig('eta_c.png')
plt.close()