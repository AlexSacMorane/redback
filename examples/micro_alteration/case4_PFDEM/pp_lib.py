# this code presents the functions available for postprocessing the simulation.
# last update: 09/09/2026
# authors: Alexandre Sac-Morane, alexandre.sac-morane@enpc.fr

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Librairies
# ------------------------------------------------------------------------------------------------------------------------------------------ #

import numpy as np

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------------------------------------------------------------ #

def compute_cog(dict_user, dict_sample):
    '''
    Compute the center of gravity of the grain
    '''
    # initialization
    cog = np.array([0, 0])
    sum_eta = 0
    # iterate on the mesh
    for i_x in range(dict_user['n_mesh_x']):
        for i_y in range(dict_user['n_mesh_y']):
            # compute the center of gravity
            cog = cog + np.array([dict_sample['x_L'][i_x], dict_sample['y_L'][i_y]])*dict_sample['eta_1_map'][-1-i_y, i_x]
            sum_eta = sum_eta + dict_sample['eta_1_map'][-1-i_y, i_x]
    # do the mean and save it
    dict_sample['L_cog_y'].append(cog[1]/sum_eta)

# ------------------------------------------------------------------------------------------------------------------------------------------ #

def compute_as_max(dict_user, dict_sample):
    '''
    Save the maximum solid activity (relative to the pressure transmitted)
    '''
    dict_sample['L_as_max'].append(np.max(dict_sample['as_map']))

# ------------------------------------------------------------------------------------------------------------------------------------------ #

def compute_distribution_mass(dict_user, dict_sample):
    '''
    Compute the distribution of the mass between the solute c and the solid eta_1
    '''
    dict_sample['L_s_c'].append(np.sum(dict_sample['c_map']))
    dict_sample['L_s_eta_1'].append(np.sum(dict_sample['eta_1_map']))
    