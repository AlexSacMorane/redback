# this code presents the functions employed in the preparation of the Phase-Field step.
# last update: 09/09/2026
# authors: Alexandre Sac-Morane, alexandre.sac-morane@enpc.fr

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Librairies
# ------------------------------------------------------------------------------------------------------------------------------------------ #

import numpy as np
from scipy.ndimage import binary_dilation

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------------------------------------------------------------ #

def update_phase_map(dict_user, dict_sample):
    '''
    Apply a rigid body motion to the phase map of the grain to consider its displacement.
    '''
    # save the previous phase map
    eta_1_map_prev = dict_sample['eta_1_map'].copy()
    # initialize the new phase map
    dict_sample['eta_1_map'] = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))

    # iterate over the y coordinates
    for i_y in range(len(dict_sample['y_L'])):
        # compute the coordinate y in the new configuration
        y = dict_sample['y_L'][i_y]
        # initialize the index for the corresponding position in the previous configuration
        i_y_old = 0

        # grain displacement in the +y direction
        if dict_sample['displacement'][1] < 0:
            if y-dict_sample['displacement'][1] <= dict_sample['y_L'][-1]:
                # look for two neareast nodes in the previous configuration
                while not (dict_sample['y_L'][i_y_old] <= y-dict_sample['displacement'][1] and y-dict_sample['displacement'][1] <= dict_sample['y_L'][i_y_old+1]):
                    i_y_old = i_y_old + 1
                # interpolate
                dict_sample['eta_1_map'][-1-i_y, :] = (eta_1_map_prev[-1-(i_y_old+1), :] - eta_1_map_prev[-1-i_y_old, :])/(dict_sample['y_L'][i_y_old+1] - dict_sample['y_L'][i_y_old])*\
                                                       (y-dict_sample['displacement'][1] - dict_sample['y_L'][i_y_old]) + eta_1_map_prev[-1-i_y_old, :]

        # grain displacement in the -y direction
        elif dict_sample['displacement'][1] > 0:
            if dict_sample['y_L'][0] <= y-dict_sample['displacement'][1]:
                # look for two neareast nodes in the previous configuration
                while not (dict_sample['y_L'][i_y_old] <= y-dict_sample['displacement'][1] and y-dict_sample['displacement'][1] <= dict_sample['y_L'][i_y_old+1]):
                    i_y_old = i_y_old + 1
                # interpolate
                dict_sample['eta_1_map'][-1-i_y, :] = (eta_1_map_prev[-1-(i_y_old+1), :] - eta_1_map_prev[-1-i_y_old, :])/(dict_sample['y_L'][i_y_old+1] - dict_sample['y_L'][i_y_old])*\
                                                        (y-dict_sample['displacement'][1] - dict_sample['y_L'][i_y_old]) + eta_1_map_prev[-1-i_y_old, :]
        else :
            dict_sample['eta_1_map'][-1-i_y, :] = eta_1_map_prev[-1-i_y, :]

# ---------------------------------------------------------------------

def update_solute_map(dict_user, dict_sample):
    '''
    Update the solute map by pushing out the solute located in the solid.
    '''
    # compute the pore space
    pore_map = np.array(np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x'])), dtype = bool)
    # iterate on the mesh
    for i_y in range(len(dict_sample['y_L'])):
        for i_x in range(len(dict_sample['x_L'])):
            if dict_sample['eta_1_map'][-1-i_y, i_x] < 0.1 and dict_sample['y_L'][i_y] > 0: # out of the grain
                pore_map[-1-i_y, i_x] = True
            elif dict_sample['eta_1_map'][-1-i_y, i_x] > 0.1 and dict_sample['y_L'][i_y] <= 0: # in the contact
                pore_map[-1-i_y, i_x] = True
            else :
                pore_map[-1-i_y, i_x] = False
    # dilate the pore space
    dilated_pore_map = binary_dilation(pore_map, np.array(np.ones((5, 5)), dtype=bool))

    # save the previous solute map
    c_map_prev = dict_sample['c_map'].copy()
    # initialize the new solute map
    dict_sample['c_map'] = dict_sample['c_map'].copy()

    # iterate on the mesh
    for i_y in range(len(dict_sample['y_L'])):
        for i_x in range(len(dict_sample['x_L'])):
            # push solute out of the solid
            if not dilated_pore_map[i_y, i_x] and c_map_prev[i_y, i_x] > 1: # threshold value
                solute_moved = False
                size_window = 1
                # compute solute to move
                solute_to_move = c_map_prev[i_y, i_x] - 1
                while not solute_moved :
                    i_window = 0
                    while not solute_moved and i_window <= size_window:
                        n_node_available = 0

                        #Look to move horizontaly and vertically
                        if i_window == 0 :
                            top_available = False
                            down_available = False
                            left_available = False
                            right_available = False
                            #to the top
                            if i_y - size_window > 0:
                                top_available = dilated_pore_map[i_y-size_window, i_x]
                                if dilated_pore_map[i_y-size_window, i_x] :
                                    n_node_available = n_node_available + 1
                            #to the down
                            if i_y + size_window < len(dict_sample['y_L']):
                                down_available = dilated_pore_map[i_y+size_window, i_x]
                                if dilated_pore_map[i_y+size_window, i_x] :
                                    n_node_available = n_node_available + 1
                            #to the left
                            if i_x - size_window > 0:
                                left_available = dilated_pore_map[i_y, i_x-size_window]
                                if dilated_pore_map[i_y, i_x-size_window] :
                                    n_node_available = n_node_available + 1
                            #to the right
                            if i_x + size_window < len(dict_sample['x_L']):
                                right_available = dilated_pore_map[i_y, i_x+size_window]
                                if dilated_pore_map[i_y, i_x+size_window] :
                                    n_node_available = n_node_available + 1

                            #move solute if et least one node is available
                            if n_node_available != 0 :
                                #to the top
                                if top_available:
                                    dict_sample['c_map'][i_y-size_window, i_x] = dict_sample['c_map'][i_y-size_window, i_x] + solute_to_move/n_node_available
                                #to the down
                                if down_available:
                                    dict_sample['c_map'][i_y+size_window, i_x] = dict_sample['c_map'][i_y+size_window, i_x] + solute_to_move/n_node_available
                                #to the left
                                if left_available:
                                    dict_sample['c_map'][i_y, i_x-size_window] = dict_sample['c_map'][i_y, i_x-size_window] + solute_to_move/n_node_available
                                #to the right
                                if right_available:
                                    dict_sample['c_map'][i_y, i_x+size_window] = dict_sample['c_map'][i_y, i_x+size_window] + solute_to_move/n_node_available
                                dict_sample['c_map'][i_y, i_x] = 1
                                solute_moved = True

                        #Look to move diagonally
                        else :
                            top_min_available = False
                            top_max_available = False
                            down_min_available = False
                            down_max_available = False
                            left_min_available = False
                            left_max_available = False
                            right_min_available = False
                            right_max_available = False
                            #to the top
                            if i_y - size_window > 0:
                                if i_x - i_window > 0 :
                                    top_min_available = dilated_pore_map[i_y-size_window, i_x-i_window]
                                    if dilated_pore_map[i_y-size_window, i_x-i_window] :
                                        n_node_available = n_node_available + 1
                                if i_x + i_window < len(dict_sample['x_L']):
                                    top_max_available = dilated_pore_map[i_y-size_window, i_x+i_window]
                                    if dilated_pore_map[i_y-size_window, i_x+i_window] :
                                        n_node_available = n_node_available + 1
                            #to the down
                            if i_y + size_window < len(dict_sample['y_L']):
                                if i_x - i_window > 0 :
                                    down_min_available = dilated_pore_map[i_y+size_window, i_x-i_window]
                                    if dilated_pore_map[i_y+size_window, i_x-i_window] :
                                        n_node_available = n_node_available + 1
                                if i_x + i_window < len(dict_sample['x_L']):
                                    down_max_available = dilated_pore_map[i_y+size_window, i_x+i_window]
                                    if dilated_pore_map[i_y+size_window, i_x+i_window] :
                                        n_node_available = n_node_available + 1
                            #to the left
                            if i_x - size_window > 0:
                                if i_y - i_window > 0 :
                                    left_min_available = dilated_pore_map[i_y-i_window, i_x-size_window]
                                    if dilated_pore_map[i_y-i_window, i_x-size_window] :
                                        n_node_available = n_node_available + 1
                                if i_y + i_window < len(dict_sample['y_L']):
                                    left_max_available = dilated_pore_map[i_y+i_window, i_x-size_window]
                                    if dilated_pore_map[i_y+i_window, i_x-size_window] :
                                        n_node_available = n_node_available + 1
                            #to the right
                            if i_x + size_window < len(dict_sample['x_L']):
                                if i_x - i_window > 0 :
                                    right_min_available = dilated_pore_map[i_y-i_window, i_x+size_window]
                                    if dilated_pore_map[i_y-i_window, i_x+size_window] :
                                        n_node_available = n_node_available + 1
                                if i_y + i_window < len(dict_sample['y_L']):
                                    right_max_available = dilated_pore_map[i_y+i_window, i_x+size_window]
                                    if dilated_pore_map[i_y+i_window, i_x+size_window] :
                                        n_node_available = n_node_available + 1

                            #move solute if et least one node is available
                            if n_node_available != 0 :
                                #to the top
                                if top_min_available:
                                    dict_sample['c_map'][i_y-size_window, i_x-i_window] = dict_sample['c_map'][i_y-size_window, i_x-i_window] + solute_to_move/n_node_available
                                if top_max_available:
                                    dict_sample['c_map'][i_y-size_window, i_x+i_window] = dict_sample['c_map'][i_y-size_window, i_x+i_window] + solute_to_move/n_node_available
                                #to the down
                                if down_min_available:
                                    dict_sample['c_map'][i_y+size_window, i_x-i_window] = dict_sample['c_map'][i_y+size_window, i_x-i_window] + solute_to_move/n_node_available
                                if down_max_available:
                                    dict_sample['c_map'][i_y+size_window, i_x+i_window] = dict_sample['c_map'][i_y+size_window, i_x+i_window] + solute_to_move/n_node_available
                                #to the left
                                if left_min_available:
                                    dict_sample['c_map'][i_y-i_window, i_x-size_window] = dict_sample['c_map'][i_y-i_window, i_x-size_window] + solute_to_move/n_node_available
                                if left_max_available:
                                    dict_sample['c_map'][i_y+i_window, i_x-size_window] = dict_sample['c_map'][i_y+i_window, i_x-size_window] + solute_to_move/n_node_available
                                #to the right
                                if right_min_available:
                                    dict_sample['c_map'][i_y-i_window, i_x+size_window] = dict_sample['c_map'][i_y-i_window, i_x+size_window] + solute_to_move/n_node_available
                                if right_max_available:
                                    dict_sample['c_map'][i_y+i_window, i_x+size_window] = dict_sample['c_map'][i_y+i_window, i_x+size_window] + solute_to_move/n_node_available
                                dict_sample['c_map'][i_y, i_x] = 1
                                solute_moved = True
                        i_window = i_window + 1
                    size_window = size_window + 1

            # push solute in of the solid
            if not dilated_pore_map[i_y, i_x] and c_map_prev[i_y, i_x] < 1: # threshold value
                solute_moved = False
                size_window = 1
                # compute solute to move
                solute_to_move = 1 - c_map_prev[i_y, i_x]
                while not solute_moved :
                    i_window = 0
                    while not solute_moved and i_window <= size_window:
                        n_node_available = 0

                        #Look to move horizontaly and vertically
                        if i_window == 0 :
                            top_available = False
                            down_available = False
                            left_available = False
                            right_available = False
                            #to the top
                            if i_y - size_window > 0:
                                top_available = dilated_pore_map[i_y-size_window, i_x]
                                if dilated_pore_map[i_y-size_window, i_x] :
                                    n_node_available = n_node_available + 1
                            #to the down
                            if i_y + size_window < len(dict_sample['y_L']):
                                down_available = dilated_pore_map[i_y+size_window, i_x]
                                if dilated_pore_map[i_y+size_window, i_x] :
                                    n_node_available = n_node_available + 1
                            #to the left
                            if i_x - size_window > 0:
                                left_available = dilated_pore_map[i_y, i_x-size_window]
                                if dilated_pore_map[i_y, i_x-size_window] :
                                    n_node_available = n_node_available + 1
                            #to the right
                            if i_x + size_window < len(dict_sample['x_L']):
                                right_available = dilated_pore_map[i_y, i_x+size_window]
                                if dilated_pore_map[i_y, i_x+size_window] :
                                    n_node_available = n_node_available + 1

                            #move solute if et least one node is available
                            if n_node_available != 0 :
                                #to the top
                                if top_available:
                                    dict_sample['c_map'][i_y-size_window, i_x] = dict_sample['c_map'][i_y-size_window, i_x] - solute_to_move/n_node_available
                                #to the down
                                if down_available:
                                    dict_sample['c_map'][i_y+size_window, i_x] = dict_sample['c_map'][i_y+size_window, i_x] - solute_to_move/n_node_available
                                #to the left
                                if left_available:
                                    dict_sample['c_map'][i_y, i_x-size_window] = dict_sample['c_map'][i_y, i_x-size_window] - solute_to_move/n_node_available
                                #to the right
                                if right_available:
                                    dict_sample['c_map'][i_y, i_x+size_window] = dict_sample['c_map'][i_y, i_x+size_window] - solute_to_move/n_node_available
                                dict_sample['c_map'][i_y, i_x] = 1
                                solute_moved = True

                        #Look to move diagonally
                        else :
                            top_min_available = False
                            top_max_available = False
                            down_min_available = False
                            down_max_available = False
                            left_min_available = False
                            left_max_available = False
                            right_min_available = False
                            right_max_available = False
                            #to the top
                            if i_y - size_window > 0:
                                if i_x - i_window > 0 :
                                    top_min_available = dilated_pore_map[i_y-size_window, i_x-i_window]
                                    if dilated_pore_map[i_y-size_window, i_x-i_window] :
                                        n_node_available = n_node_available + 1
                                if i_x + i_window < len(dict_sample['x_L']):
                                    top_max_available = dilated_pore_map[i_y-size_window, i_x+i_window]
                                    if dilated_pore_map[i_y-size_window, i_x+i_window] :
                                        n_node_available = n_node_available + 1
                            #to the down
                            if i_y + size_window < len(dict_sample['y_L']):
                                if i_x - i_window > 0 :
                                    down_min_available = dilated_pore_map[i_y+size_window, i_x-i_window]
                                    if dilated_pore_map[i_y+size_window, i_x-i_window] :
                                        n_node_available = n_node_available + 1
                                if i_x + i_window < len(dict_sample['x_L']):
                                    down_max_available = dilated_pore_map[i_y+size_window, i_x+i_window]
                                    if dilated_pore_map[i_y+size_window, i_x+i_window] :
                                        n_node_available = n_node_available + 1
                            #to the left
                            if i_x - size_window > 0:
                                if i_y - i_window > 0 :
                                    left_min_available = dilated_pore_map[i_y-i_window, i_x-size_window]
                                    if dilated_pore_map[i_y-i_window, i_x-size_window] :
                                        n_node_available = n_node_available + 1
                                if i_y + i_window < len(dict_sample['y_L']):
                                    left_max_available = dilated_pore_map[i_y+i_window, i_x-size_window]
                                    if dilated_pore_map[i_y+i_window, i_x-size_window] :
                                        n_node_available = n_node_available + 1
                            #to the right
                            if i_x + size_window < len(dict_sample['x_L']):
                                if i_x - i_window > 0 :
                                    right_min_available = dilated_pore_map[i_y-i_window, i_x+size_window]
                                    if dilated_pore_map[i_y-i_window, i_x+size_window] :
                                        n_node_available = n_node_available + 1
                                if i_y + i_window < len(dict_sample['y_L']):
                                    right_max_available = dilated_pore_map[i_y+i_window, i_x+size_window]
                                    if dilated_pore_map[i_y+i_window, i_x+size_window] :
                                        n_node_available = n_node_available + 1

                            #move solute if et least one node is available
                            if n_node_available != 0 :
                                #to the top
                                if top_min_available:
                                    dict_sample['c_map'][i_y-size_window, i_x-i_window] = dict_sample['c_map'][i_y-size_window, i_x-i_window] - solute_to_move/n_node_available
                                if top_max_available:
                                    dict_sample['c_map'][i_y-size_window, i_x+i_window] = dict_sample['c_map'][i_y-size_window, i_x+i_window] - solute_to_move/n_node_available
                                #to the down
                                if down_min_available:
                                    dict_sample['c_map'][i_y+size_window, i_x-i_window] = dict_sample['c_map'][i_y+size_window, i_x-i_window] - solute_to_move/n_node_available
                                if down_max_available:
                                    dict_sample['c_map'][i_y+size_window, i_x+i_window] = dict_sample['c_map'][i_y+size_window, i_x+i_window] - solute_to_move/n_node_available
                                #to the left
                                if left_min_available:
                                    dict_sample['c_map'][i_y-i_window, i_x-size_window] = dict_sample['c_map'][i_y-i_window, i_x-size_window] - solute_to_move/n_node_available
                                if left_max_available:
                                    dict_sample['c_map'][i_y+i_window, i_x-size_window] = dict_sample['c_map'][i_y+i_window, i_x-size_window] - solute_to_move/n_node_available
                                #to the right
                                if right_min_available:
                                    dict_sample['c_map'][i_y-i_window, i_x+size_window] = dict_sample['c_map'][i_y-i_window, i_x+size_window] - solute_to_move/n_node_available
                                if right_max_available:
                                    dict_sample['c_map'][i_y+i_window, i_x+size_window] = dict_sample['c_map'][i_y+i_window, i_x+size_window] - solute_to_move/n_node_available
                                dict_sample['c_map'][i_y, i_x] = 1
                                solute_moved = True
                        i_window = i_window + 1
                    size_window = size_window + 1

# ---------------------------------------------------------------------

def characterize_contact(dict_user, dict_sample):
    '''
    Characterize the contact between the grain and the wall.
    '''
    # recompute the contact, start the investigation from the bottom
    i_y = 0
    # compute the contact volume/size at this current coordinate
    contact_volume = dict_sample['eta_1_map'][-1-i_y,:].sum()*(dict_sample['x_L'][1]-dict_sample['x_L'][0])*(dict_sample['y_L'][1]-dict_sample['y_L'][0])*1
    contact_size = dict_sample['eta_1_map'][-1-i_y,:].sum()*(dict_sample['x_L'][1]-dict_sample['x_L'][0])
    # save the previous contact volume/size
    contact_volume_prev = 0
    contact_size_prev = 0
    # compare the contact volume to the imposed solicitation
    while contact_volume < dict_user['imposed_contact_volume']:
        # iterate over the y coordinates
        i_y = i_y + 1
        # update the contact volume and the previous contact volume
        contact_volume_prev = contact_volume
        contact_volume = contact_volume + dict_sample['eta_1_map'][-1-i_y,:].sum()*(dict_sample['x_L'][1]-dict_sample['x_L'][0])*(dict_sample['y_L'][1]-dict_sample['y_L'][0])*1
        # update the contact size and the previous contact size
        contact_size_prev = contact_size
        contact_size = dict_sample['eta_1_map'][-1-i_y,:].sum()*(dict_sample['x_L'][1]-dict_sample['x_L'][0])
    # interpolate the contact surface 
    dict_sample['contact_surface'] = (contact_size_prev + (dict_user['imposed_contact_volume']-contact_volume_prev)/(contact_volume-contact_volume_prev)*(contact_size-contact_size_prev))*1
    # save the initial contact surface
    if 'initial_contact_surface' not in dict_sample.keys():
        dict_sample['initial_contact_surface'] = dict_sample['contact_surface']    
    
# ---------------------------------------------------------------------

def compute_as_map(dict_user, dict_sample):
    '''
    Compute the solid activity map.
    '''
    # initialize the solid activity map
    dict_sample['as_map'] = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))
    # iterate on the mesh
    for i_y in range(len(dict_sample['y_L'])):
        for i_x in range(len(dict_sample['x_L'])):
            if dict_sample['eta_1_map'][-1-i_y, i_x] > 0.1 and dict_sample['y_L'][i_y] <= 0: # in the contact
                dict_sample['as_map'][-1-i_y, i_x] = 1.2**(dict_sample['initial_contact_surface']/dict_sample['contact_surface'])
            else : # not in the contact
                dict_sample['as_map'][-1-i_y, i_x] = 1

# ---------------------------------------------------------------------

def compute_diffusivity_map(dict_user, dict_sample):
    '''
    Compute the diffusivity map.

    It is assumed that the solute diffusion coefficient is larger in the pore than in the contact space.
    '''      
    # compute
    contact_map = np.array(np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x'])), dtype = bool)
    pore_map =  np.array(np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x'])), dtype = bool)
    # iterate over x and y
    for i_y in range(len(dict_sample['y_L'])):
        for i_x in range(len(dict_sample['x_L'])):
            if dict_sample['eta_1_map'][-1-i_y, i_x] < dict_user['eta_contact_box_detection'] and dict_sample['y_L'][i_y] > 0: # out of the grain
                contact_map[-1-i_y, i_x] = True
                pore_map[-1-i_y, i_x] = True
            elif dict_sample['eta_1_map'][-1-i_y, i_x] > dict_user['eta_contact_box_detection'] and dict_sample['y_L'][i_y] <= 0: # in the contact
                contact_map[-1-i_y, i_x] = True
                pore_map[-1-i_y, i_x] = False
            else :
                contact_map[-1-i_y, i_x] = False
                pore_map[-1-i_y, i_x] = False

    # dilate the contact map
    dilated_contact_map = binary_dilation(contact_map, np.array(np.ones((5, 5)), dtype=bool))

    #compute the map of the solute diffusion coefficient
    dict_sample['kc_map'] = dict_user['D_solute']*dilated_contact_map + 99*dict_user['D_solute']*pore_map