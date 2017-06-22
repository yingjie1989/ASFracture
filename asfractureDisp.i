# Wave propogation in 1D using Newmark time integration
#
# The test is for an 1D bar element of length 4m  fixed on one end
# with a sinusoidal pulse dirichlet boundary condition applied to the other end.
# beta and gamma are Newmark time integration parameters
# The equation of motion in terms of matrices is:
#
# M*accel +  K*disp = 0
#
# Here M is the mass matrix, K is the stiffness matrix
#
# This equation is equivalent to:
#
# density*accel + Div Stress= 0
#
# The first term on the left is evaluated using the Inertial force kernel
# The last term on the left is evaluated using StressDivergenceTensors
#
# The displacement at the second, third and fourth node at t = 0.1 are
# -8.021501116638234119e-02, 2.073994362053969628e-02 and  -5.045094181261772920e-03, respectively

#[Mesh]
#  type = GeneratedMesh
#  dim = 2
#  nx = 1
#  ny = 4
#  nz = 1
#  xmin = 0.0
#  xmax = 0.1
#  ymin = 0.0
#  ymax = 4.0
#  zmin = 0.0
#  zmax = 0.1
#[]

[Mesh]
  file = /home/yl306/ASFracture/mesh/rock.msh
  parallel_type = DISTRIBUTED
[]


[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[AuxVariables]
  [./vel_x]
  [../]
  [./accel_x]
  [../]
  [./vel_y]
  [../]
  [./accel_y]
  [../]
  [./vel_z]
  [../]
  [./accel_z]
  [../]
  [./pre_wave]
  [../]
#  [./stress_xx]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./strain_xx]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#   [./stress_yy]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./strain_yy]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./stress_xy]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./strain_xy]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
  [./stress_h]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_max]
    order = CONSTANT
    family = MONOMIAL
  [../]
#  [./stress_min]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./stress_VonMises]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
  [./c_in]
  [../]
[]

[Kernels]
# Dynamics Equation
  [./solid_x]
     type = StressDivergencePFFracTensors
     variable = disp_x
     displacements = 'disp_x disp_y disp_z'
     component = 0
     c = c_in
  [../]
  [./solid_y]
     type = StressDivergencePFFracTensors
     variable = disp_y
     displacements = 'disp_x disp_y disp_z'
     component = 1
     c = c_in
  [../]
  [./solid_z]
     type = StressDivergencePFFracTensors
     variable = disp_z
     displacements = 'disp_x disp_y disp_z'
     component = 2
     c = c_in
  [../]
  [./inertia_x]
    type = InertialForce
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = 0.3025
    gamma = 0.6
    eta=0.0
  [../]
  [./inertia_y]
    type = InertialForce
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = 0.3025
    gamma = 0.6
    eta=0.0
  [../]
  [./inertia_z]
    type = InertialForce
    variable = disp_z
    velocity = vel_z
    acceleration = accel_z
    beta = 0.3025
    gamma = 0.6
    eta = 0.0
  [../]
[]

[AuxKernels]
  [./accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = 0.3025
    execute_on = timestep_end
  [../]
  [./vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = 0.6
    execute_on = timestep_end
  [../]
  [./accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = 0.3025
    execute_on = timestep_end
  [../]
  [./vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = 0.6
    execute_on = timestep_end
  [../]
  [./accel_z]
    type = NewmarkAccelAux
    variable = accel_z
    displacement = disp_z
    velocity = vel_z
    beta = 0.3025
    execute_on = timestep_end
  [../]
  [./vel_z]
    type = NewmarkVelAux
    variable = vel_z
    acceleration = accel_z
    gamma = 0.6
    execute_on = timestep_end
  [../]
#  [./stress_xx]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    variable = stress_xx
#    index_i = 0
#    index_j = 0
#  [../]
#  [./strain_xx]
#    type = RankTwoAux
#    rank_two_tensor = total_strain
#    variable = strain_xx
#    index_i = 0
#    index_j = 0
#  [../]
#   [./stress_yy]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    variable = stress_yy
#    index_i = 1
#    index_j = 1
#  [../]
#  [./strain_yy]
#    type = RankTwoAux
#    rank_two_tensor = total_strain
#    variable = strain_yy
#    index_i = 1
#    index_j = 1
#  [../]
#  [./stress_xy]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    variable = stress_xy
#    index_i = 0
#    index_j = 1
#  [../]
#  [./strain_xy]
#    type = RankTwoAux
#    rank_two_tensor = total_strain
#    variable = strain_xy
#    index_i = 0
#    index_j = 1
#  [../]
  [./stress_h]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = stress_h
    scalar_type = Hydrostatic
  [../]
  [./stress_max]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = stress_max
    scalar_type = MaxPrincipal
  [../]
#  [./stress_min]
#    type = RankTwoScalarAux
#    rank_two_tensor = stress
#    variable = stress_min
#    scalar_type = MinPrincipal
#  [../]
#  [./stress_VonMises]
#    type = RankTwoScalarAux
#    rank_two_tensor = stress
#    variable = stress_VonMises
#    scalar_type = VonMisesStress
#  [../] 
[]


[BCs]
  [./top_z]
    type = CoupledNeumannBC
    variable = disp_z
    boundary = top2
    some_var = pre_wave
    Reyold = -1.0
  [../]
#  [./top_x]
#    type = DirichletBC
#    variable = disp_x
#    boundary = 4
#    value=0.0
#  [../]
#  [./top_z]
#    type = DirichletBC
#    variable = disp_z
#    boundary = top
#    value=0.0
#  [../]
  [./inner_x]
    type = CoupledNeumannVectorBC
    variable = disp_x
    boundary = inner
    some_var_x = pre_wave
    Reyold = -1.0
  [../]
#  [./right_z]
#    type = DirichletBC
#    variable = disp_z
#    boundary = 3
#    value=0.0
#  [../]
  [./inner_y]
    type = CoupledNeumannVectorBC
    variable = disp_y
    boundary = inner
    some_var_y = pre_wave
    Reyold = -1.0
  [../]
#  [./left_z]
#    type = DirichletBC
#    variable = disp_z
#    boundary = left
#    value=0.0
#  [../]
#  [./front_x]
#    type = DirichletBC
#    variable = disp_x
#    boundary = front
#    value=0.0
#  [../]
#  [./front_z]
#    type = DirichletBC
#    variable = disp_z
#    boundary = front
#    value=0.0
#  [../]
#  [./back_x]
#    type = DirichletBC
#    variable = disp_x
#    boundary = back
#    value=0.0
#  [../]
#  [./back_z]
#    type = DirichletBC
#    variable = disp_z
#    boundary = back
#    value=0.0
#  [../]
#  [./bottom_x]
#    type = DirichletBC
#    variable = disp_x
#    boundary = 2
#    value=0.0
#  [../]
#  [./bottom_z]
#    type = DirichletBC
#    variable = disp_z
#    boundary = bottom
#    value=0.0
#  [../]
  [./bottom_z]
    type = CoupledNeumannBC
    variable = disp_z
    boundary = bottom2
    some_var = pre_wave
    Reyold = 1.0
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    fill_method = symmetric_isotropic
    C_ijkl = '0.01305 0.01073'
  [../]

  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y disp_z'
  [../]

  [./elastic]
       type = LinearIsoElasticPFDamage
       c = c_in
       kdamage = 1e-6
  [../]
  [./density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '1.995e-3'
  [../]
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]
[]


[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      31                  preonly       lu           1'

  start_time = 0
  end_time = 3.00
  l_max_its = 200
#  nl_max_its = 10
  
  l_tol = 1e-12
  nl_rel_tol = 1e-12
  dt = 5e-3
[]

[Outputs]
  csv = true
  [./console]
        type = Console
        max_rows = 10
  [../]
  [./exodus]
        type = Exodus
        interval = 5
        file_base = DispOut
  [../]
[]

[MultiApps]
  [./sub_damage]
    type = TransientMultiApp
    app_type = ASFracture
    execute_on = timestep_begin
    positions = '0.0 0.0 0.0'
    input_files = asfractureDamage.i
  [../]
[]

[Transfers]
  [./from_sub_c]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub_damage
    source_variable = c
    variable = c_in
    execute_on = 'timestep_begin'
  [../]
  [./to_sub_x]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub_damage
    source_variable = disp_x
    variable = disp_x_in
    execute_on = 'timestep_end'
  [../]
  [./to_sub_y]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub_damage
    source_variable = disp_y
    variable = disp_y_in
    execute_on = 'timestep_end'
  [../]
  [./to_sub_z]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub_damage
    source_variable = disp_z
    variable = disp_z_in
    execute_on = 'timestep_end'
  [../]
[]

