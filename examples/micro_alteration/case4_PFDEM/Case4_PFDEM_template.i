# this code is a part of the PFDEM approache
# the chemical destabilization (dissolution/precipitation) is estimated
# the "Phase Field" module of MOOSE should be used for the resolution
# last update: 09/09/2026
# authors: Alexandre Sac-Morane, alexandre.sac-morane@enpc.fr

# Domain and mesh descriptions
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 
  ny = 
  xmin = 
  xmax = 
  ymin = 
  ymax = 
  elem_type = QUAD4
[]

# Definition of the variables
[Variables]
  # phase variable
  [./eta_1]
    order = FIRST
    family = LAGRANGE
    outputs = exodus
    [./InitialCondition]
      type = FunctionIC
      function = eta_1_txt
    [../]
  [../]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = c_txt
    [../]
  [../]
[]

# definition of the problem
[Kernels]
  #
  # variable eta_1
  #
  # time derivative
  [./deta_1dt]
    type = TimeDerivative
    variable = eta_1
  [../]
  # free energy and tilting energy
  [./ACBulk_eta_1]
    type = AllenCahn
    variable = eta_1
    mob_name = L_eta
    f_name = g_eta_1
  [../]
  # diffusive term
  [./ACInterface_eta]
    type = ACInterface
    variable = eta_1
    mob_name = L_eta
    kappa_name = kappa_eta
  [../]

  #
  # variable c
  #
  [./dcdt]
    type = TimeDerivative
    variable = c
  [../]
  [./eta_1_c]
    type = CoefCoupledTimeDerivative
    v = 'eta_1'
    variable = c
    coef = 1
  [../]
  [./c_diffusion]
    type = ACInterface
    kappa_name = kc
    mob_name = L_c
    variable = c
  [../]
[]

# Material description
[Materials]
  # constants of the problem
  [./consts]
    # L_eta can be changed to modify the reaction kinetics
    # The diffusion coefficient (kappa_eta) induce also a dissolution of the phase.
    # the mechanisms due to the tilting energy should be larger than the mechanisms induced by the diffusion
    type = GenericConstantMaterial
    prop_names  = 'L_eta kappa_eta L_c'
    prop_values = 
  [../]
  # map of the diffusivity in the domain
  [./kcmap]
    type = GenericFunctionMaterial
    prop_names = kc
    prop_values = kc_txt
    outputs = exodus
  [../]
  # map of the solid acitvity in the domain
  [./asmap]
    type = GenericFunctionMaterial
    prop_names = as
    prop_values = as_txt
    outputs = exodus
  [../]
  # energy of the phase variable
  # this energy is the free energy + the tilting energy
  [./free_energy_and_e]
    type = DerivativeParsedMaterial
    property_name = g_eta_1
    coupled_variables = 'eta_1 c'
    material_property_names = 'F(eta_1) E(eta_1,c)'
    expression = 'F+E'
    enable_jit = true
    derivative_order = 2
    outputs = exodus
  [../]
  # free energy
  [./free_energy]
    type = DerivativeParsedMaterial
    property_name = F
    coupled_variables = 'eta_1'
    constant_names = 'W'
    constant_expressions = '1' # constants can be changed
    expression = 'W*(eta_1^2)*((1-eta_1)^2)'
    enable_jit = true
    derivative_order = 2
    #outputs = exodus
  [../]
  # the tilting energy
  [./Ed]
    type = DerivativeParsedMaterial
    property_name = E
    coupled_variables = 'eta_1 c'
    material_property_names = 'as'
    constant_names = 'c_eq k_diss k_prec' # the dissolution/precipication kinetics can differ
    constant_expressions =
    expression = 'if(c<c_eq*as,k_diss*as*(1-c/(c_eq*as))*(3*eta_1^2-2*eta_1^3),k_prec*as*(1-c/(c_eq*as))*(3*eta_1^2-2*eta_1^3))'
    enable_jit = true
    derivative_order = 2
    outputs = exodus
  [../]
[]

[Functions]
  [eta_1_txt]
    type = PiecewiseMultilinear
    data_file = MOOSE_simulation/eta_1_map.txt
  []
  [c_txt]
    type = PiecewiseMultilinear
    data_file = MOOSE_simulation/c_map.txt
  []
	[as_txt]
		type = PiecewiseMultilinear
		data_file = MOOSE_simulation/as_map.txt
	[]
	[kc_txt]
		type = PiecewiseMultilinear
		data_file = MOOSE_simulation/kc_map.txt
	[]
[]

[Preconditioning]
  # This preconditioner makes sure the Jacobian Matrix is fully populated. Our
  # kernels compute all Jacobian matrix entries.
  # This allows us to use the Newton solver below.
  [./SMP]
    type = SMP
    full = true
  [../]
[]

# define the resolution scheme
[Executioner]
  type = Transient
  scheme = 'bdf2'
  # Automatic differentiation provides a _full_ Jacobian in this example
  # so we can safely use NEWTON for a fast solve
  solve_type = 'NEWTON'
  l_max_its = 20
  l_tol = 1e-6
  l_abs_tol = 1e-6
  nl_max_its = 10
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  # time domain of the simulation
  start_time = 0.0
  end_time = 
  # use an adaptative time increment
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 
  [../]
[]

# Definition of the outputs
[Outputs]
  execute_on = 'initial timestep_end'
  exodus = true # general output readable in Paraview
  [console] # user output to follow the simulation
    type = Console
    execute_on = 'timestep_end'
    all_variable_norms = true
    max_rows = 5
  []
  # a vtk output is required for the parser between the methods
  [./other]
    type = VTK
    execute_on = 'initial final'
  [../]
[]
