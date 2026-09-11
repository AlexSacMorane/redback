# this code presents the functions employed during the call of the Phase-Field simulation with MOOSE.
# last update: 09/09/2026
# authors: Alexandre Sac-Morane, alexandre.sac-morane@enpc.fr

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Librairies
# ------------------------------------------------------------------------------------------------------------------------------------------ #

import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------------------------------------------------------------ #

def write_map_txt(dict_user, dict_sample, var_name):
    '''
    Write the map of a variable to a .txt file (format required for MOOSE simulation).
    '''
    # open the file
    file_to_write = open('MOOSE_simulation/' + var_name + '.txt','w')
    # x
    file_to_write.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # y
    file_to_write.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # data
    file_to_write.write('DATA\n')
    for j in range(len(dict_sample['y_L'])):
        for i in range(len(dict_sample['x_L'])):
            # variable
            file_to_write.write(str(dict_sample[var_name][-1-j,i])+'\n')
    # close the file
    file_to_write.close()

# ---------------------------------------------------------------------

def write_moose_input_file(dict_user, dict_sample):
    '''
    Write the MOOSE input file for the Phase-Field simulation.

    The file is generated from a template nammed Case4_PFDEM_template.i
    '''
    file_to_write = open('Case4_PFDEM.i','w')
    file_to_read = open('Case4_PFDEM_template.i','r')
    lines = file_to_read.readlines()
    file_to_read.close()

    # read the template and fill it if required
    j = 0
    for line in lines :
        j = j + 1
        if j == 11:
            line = line[:-1] + ' ' + str(len(dict_sample['x_L'])-1)+'\n'
        elif j == 12:
            line = line[:-1] + ' ' + str(len(dict_sample['y_L'])-1)+'\n'
        elif j == 13:
            line = line[:-1] + ' ' + str(min(dict_sample['x_L']))+'\n'
        elif j == 14:
            line = line[:-1] + ' ' + str(max(dict_sample['x_L']))+'\n'
        elif j == 15:
            line = line[:-1] + ' ' + str(min(dict_sample['y_L']))+'\n'
        elif j == 16:
            line = line[:-1] + ' ' + str(max(dict_sample['y_L']))+'\n'
        elif j == 95:
            line = line[:-1] + " '1 " + str(dict_user['kappa_eta']) + " 1'\n"
        elif j == 142:
            line = line[:-1] + " '1 " + str(dict_user['k_reac']) + ' ' + str(2*dict_user['k_reac']) + "'\n"
        elif j == 194:
            line = line[:-1] + ' ' + str(dict_user['time_PF']) + '\n'
        elif j == 198:
            line = line[:-1] + ' ' + str(dict_user['dt_PF']) + '\n'

        file_to_write.write(line)

    file_to_write.close()

# ---------------------------------------------------------------------

def adapt_pvtu(filename, index_str):
    '''
    Adapt the pvtu files
    '''
    file_to_write = open('vtk/'+filename[:-8]+index_str+'.pvtu','w')
    file_to_read = open(filename,'r')
    lines = file_to_read.readlines()
    file_to_read.close()

    # read the template and fill it if required
    for line in lines :
        if '<Piece' in line: 
            line = line[:37] + index_str + line[40:]
        file_to_write.write(line)
    file_to_write.close()

# ---------------------------------------------------------------------

def read_vtk(dict_user, dict_sample, index_str):
    '''
    Read the last vtk files to obtain data from MOOSE.
    '''
    if not dict_sample['Map_known']:
        L_XY = []
        L_L_i_XY = []

    # iterate on the proccessors used
    for i_proc in range(dict_user['n_proc']):
        # name of the file to load
        namefile = 'vtk/Case4_PFDEM_other_'+index_str+'_'+str(i_proc)+'.vtu'

        # load a vtk file as input
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(namefile)
        reader.Update()

        # Grab a scalar from the vtk file
        nodes_vtk_array = reader.GetOutput().GetPoints().GetData()
        eta_1_vtk_array = reader.GetOutput().GetPointData().GetArray("eta_1")
        c_vtk_array = reader.GetOutput().GetPointData().GetArray("c")

        #Get the coordinates of the nodes and the scalar values
        nodes_array = vtk_to_numpy(nodes_vtk_array)
        eta_1_array = vtk_to_numpy(eta_1_vtk_array)
        c_array = vtk_to_numpy(c_vtk_array)

        # map is not know
        if not dict_sample['Map_known']:
            # save the map
            L_i_XY = []
            # Must detect common zones between processors
            for i_XY in range(len(nodes_array)) :
                XY = nodes_array[i_XY]
                # Do not consider twice a point
                if list(XY) not in L_XY :
                    # search node in the mesh
                    L_search = list(abs(np.array(dict_sample['x_L']-list(XY)[0])))
                    i_x = L_search.index(min(L_search))
                    L_search = list(abs(np.array(dict_sample['y_L']-list(XY)[1])))
                    i_y = L_search.index(min(L_search))
                    # save map
                    L_XY.append(list(XY))
                    L_i_XY.append([i_x, i_y])
                    # rewrite map
                    dict_sample['eta_1_map'][-1-i_y, i_x] = eta_1_array[i_XY]
                    dict_sample['c_map'][-1-i_y, i_x] = c_array[i_XY]
                else :
                    L_i_XY.append(None)
            # Here the algorithm can be help as the mapping is known
            L_L_i_XY.append(L_i_XY)

        # map is known
        else :
            # iterate on data
            for i_XY in range(len(nodes_array)) :
                # read
                if dict_sample['L_L_i_XY_used'][i_proc][i_XY] != None :
                    i_x = dict_sample['L_L_i_XY_used'][i_proc][i_XY][0]
                    i_y = dict_sample['L_L_i_XY_used'][i_proc][i_XY][1]
                    # rewrite map
                    dict_sample['eta_1_map'][-1-i_y, i_x] = eta_1_array[i_XY]
                    dict_sample['c_map'][-1-i_y, i_x] = c_array[i_XY]
    
    if not dict_sample['Map_known']:
        # the map is known
        dict_sample['Map_known'] = True
        dict_sample['L_L_i_XY_used'] = L_L_i_XY