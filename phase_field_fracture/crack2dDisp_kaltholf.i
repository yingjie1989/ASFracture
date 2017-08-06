#Fix y displacement on the left side
[Mesh]
  type = FileMesh
  file = /home/yl306/ASFracture/mesh/kaltholf.msh
  uniform_refine = 0
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[AuxVariables]
#  [./resid_x]
#  [../]
#  [./resid_y]
#  [../]
[./vel_x]
  type = SideInitialCondition
  boundary = rightdown
[../]
[./accel_x]
[../]
[./vel_y]
[../]
[./accel_y]
[../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./d]
  [../]
[]

[Functions]
  [./tfunc]
    type = ParsedFunction
    value = t
  [../]
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = 0.25
    gamma = 0.5
    eta=0.0
  [../]
  [./inertia_y]
    type = InertialForce
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = 0.25
    gamma = 0.5
    eta=0.0
  [../]
  [./solid_x]
     type = StressDivergencePFFracTensors
     variable = disp_x
     displacements = 'disp_x disp_y'
     component = 0
     c = d
  [../]
  [./solid_y]
     type = StressDivergencePFFracTensors
     variable = disp_y
     displacements = 'disp_x disp_y'
     component = 1
     c = d
  [../]
[]

[AuxKernels]
  [./accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = 0.25
    execute_on = timestep_end
  [../]
  [./vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = 0.5
    execute_on = timestep_end
  [../]
  [./accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = 0.25
    execute_on = timestep_end
  [../]
  [./vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = 0.5
    execute_on = timestep_end
  [../]
  [./stress_xy]
    type = RankTwoAux
    variable = stress_xy
    rank_two_tensor = stress
    index_j = 0
    index_i = 1
    execute_on = timestep_end
  [../]
[]

[BCs]
  #[./xdisp]
  #  type = FunctionPresetBC
  #  variable = disp_x
  #  boundary = 2
  #  function = tfunc
  #[../]
  [./yfix]
    type = PresetBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  [../]
  #[./xfix]
  #  type = PresetBC
  #  variable = disp_x
  #  boundary = 1
  #  value = 0
  #[../]
[]

[Materials]
  [./pfbulkmat]
    type = PFFracBulkRateMaterial
    gc = 2.2e+4
  [../]

  [./elastic]
       type = CohesiveLinearIsoElasticPFDamage
       c = d
       kdamage = 0.0
       store_stress_old = true
       gc_prop_var = 'gc_prop'
       Emod = 'E'
       sigmac = 'sc'
       l = 0.0025
       p = 3
  [../]

  [./constant]
    type = GenericConstantMaterial
    prop_names = 'E sc density'
    prop_values = '1.9e+11 1.07e+9 8000'
  [../]

  [./elasticity_tensor]
    type = ComputeElasticityTensor
    young = 1.9e+11
    poisson = 0.3
#    C_ijkl = '121.0 81.0'
    fill_method = symmetric_isotropic
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
  [../]
[]

#[Postprocessors]
#  [./resid_x]
#    type = NodalSum
#    variable = resid_x
#    boundary = 2
#  [../]
#  [./resid_y]
#    type = NodalSum
#    variable = resid_y
#    boundary = 2
#  [../]
#[]

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
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      31                  preonly       lu           1'

#  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
#  petsc_options_value = 'lu superlu_dist'

  nl_rel_tol = 1e-8
  l_max_its = 30
  nl_max_its = 30

  dt = 1e-5
  dtmin = 1e-10
  start_time = 0.0
  end_time = 1
  #num_steps = 2
[]

[Outputs]
  file_base = KaltholfDisp
  interval = 10
  exodus = true
  csv = true
  gnuplot = true
[]

[MultiApps]
  [./sub_damage]
    type = TransientMultiApp
    app_type = ASFracture
    execute_on = timestep_begin
    positions = '0.0 0.0 0.0'
    input_files = crack2dDamage_kaltholf.i
  [../]
[]

[Transfers]
  [./from_sub_c]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub_damage
    source_variable = d
    variable = d
    execute_on = 'timestep_begin'
  [../]
  [./to_sub_x]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_damage
    source_variable = disp_x
    variable = disp_x_in
    execute_on = 'timestep_end'
  [../]
  [./to_sub_y]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub_damage
    source_variable = disp_y
    variable = disp_y_in
    execute_on = 'timestep_end'
  [../]
[]
