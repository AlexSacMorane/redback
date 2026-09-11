# this file is used to inform the main code about the parameters of the simulation.
# last update: 09/09/2026
# authors: Alexandre Sac-Morane, alexandre.sac-morane@enpc.fr

#------------------------------------------------------------------------------------------------------------------------------------------ #
# Librairies
#------------------------------------------------------------------------------------------------------------------------------------------#

import numpy as np

#------------------------------------------------------------------------------------------------------------------------------------------ #
# Parameters
#------------------------------------------------------------------------------------------------------------------------------------------#

def get_parameters():
    '''
    Define the parameters used in the simulation.
    '''
    # initialize the dictionaries
    dict_user = {}
    dict_sample = {}

    # ---------------------------------------------------------------------
    # Grain description
    
    # radius
    dict_user['radius'] = 1.0

    # ---------------------------------------------------------------------
    # Mesh description
  
    dict_user['x_min'] = -1.5*dict_user['radius']
    dict_user['x_max'] =  1.5*dict_user['radius']
    dict_user['y_min'] = -0.7*dict_user['radius']
    dict_user['y_max'] =  2.3*dict_user['radius']
    dict_user['n_mesh_x'] = 100
    dict_user['n_mesh_y'] = 100

    # compute the mean size of the mesh
    m_size_mesh = ((dict_user['x_max']-dict_user['x_min'])/(dict_user['n_mesh_x']-1)+\
                   (dict_user['y_max']-dict_user['y_min'])/(dict_user['n_mesh_y']-1))/2

    # compute the coordinate of the mesh
    dict_sample['x_L'] = np.linspace(dict_user['x_min'], dict_user['x_max'], dict_user['n_mesh_x'])
    dict_sample['y_L'] = np.linspace(dict_user['y_min'], dict_user['y_max'], dict_user['n_mesh_y'])
    
    # initialize the flag concerning the return map
    dict_sample['Map_known'] = False

    # ---------------------------------------------------------------------
    # Phase-Field description

    # number of mesh in the interface
    n_int = 6
    # the interface thickness
    dict_user['w_int'] = m_size_mesh*n_int

    # ---------------------------------------------------------------------
    # Solute diffusivity description

    dict_user['D_solute'] = 0.05

    # ---------------------------------------------------------------------
    # External solicitation description

    # the external force is constant.
    # the overlaping volume (which is proportional to the external force) is constant
    dict_user['imposed_contact_volume'] = 0.15

    # ---------------------------------------------------------------------
    # PFDEM description

    # set a maximum number of iterations
    dict_user['n_PFDEM_ite'] = 30
    # initialize the iteration counter
    dict_sample['i_PFDEM_ite'] = 0

    # number of processor
    dict_user['n_proc'] = 8 

    # ---------------------------------------------------------------------
    # PF description

    # threshold value of eta for the contact
    dict_user['eta_contact_box_detection'] = 0.1

    # the gradient coefficient for the phase variable
    dict_user['kappa_eta'] = 1*(m_size_mesh*6)**2/10

    # chemical reaction kinetics (similar for dissolution and precipitation)
    dict_user['k_reac'] = 1*(0.01)/(m_size_mesh) # ed_j = ed_i*m_i/m_j

    # duration of the PF iterations (should not induce a big change of the shape to prevent instabilities)
    dict_user['time_PF'] = 4
    # time increment of the PF simulation
    dict_user['dt_PF'] = 0.1

    # ---------------------------------------------------------------------
    # trackers initialization

    dict_sample['L_cog_y'] = []
    dict_sample['L_as_max'] = []
    dict_sample['L_s_c'] = []
    dict_sample['L_s_eta_1'] = []

    return dict_user, dict_sample