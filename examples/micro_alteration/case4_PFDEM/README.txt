# last update: 09/09/2026
# authors: Alexandre Sac-Morane, alexandre.sac-morane@enpc.fr

The global approach is based on a first-layer code that is able to generate, call, and read second-layer codes that are used for the Phase-Field and the Discrete Element Model.

This first-layer code is the python script 'main.py'. The simulation can be launched by calling this script in the terminal.
The input parameters of the simulation are available in the python script 'parameters.py'.
The scripts 'x_lib.py' are different functions required in the step 'x' of the PFDEM algorithm.
The file 'Case4_PFDEM_template.i' is a Template used by the script 'main.py' to generate the input required for the Phase-Field simulation. 

More details on the algorithm are available online.