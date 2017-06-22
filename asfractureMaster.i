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
  file = /home/yl306/ASFracture/mesh/fluid.msh
#  block_id = '0'
  parallel_type = DISTRIBUTED
#  displacements = 'disp_x disp_y'
[]


[Variables]
  [./p]
  [../]
[]

[AuxVariables]
  [./vel_p]
  [../]
  [./accel_p]
  [../]
  [./accel_x_in]
  [../]
  [./accel_y_in]
  [../]
  [./accel_z_in]
  [../]
[]

[Kernels]
  [./inertia_p]
    type = InertialForce
    variable = p
    velocity = vel_p
    acceleration = accel_p
    beta = 0.3025
    gamma = 0.6
    eta=0.0
  [../]

  [./diff_p]  
    type = Diffusion_D
    variable = p
    Diffusivity = Diff
  [../]

  [./source_p]
     type = SourceMonopole
     variable = p
     coord = '0.0 0.0 5.021'
     size  = 0.01
     fL = 8.33e-2
     t1 = 0.07
     tRT = 0.01
     tL  = 0.8
     p0  = 0.05
     d1  = 9
     upcoeff = 12.2189
     downcoeff = 0.9404
     rho_c = 1e-3
  [../]  
[]

#[DiracKernels]
#  [./monopole_source]
#    type = MonopoleDirac
#     variable = p
#     point = '0.0 0.0 5.02'
#     dim  = 3
#     fL = 8.33e-2
#     t1 = 0.07
#     tRT = 0.01
#     tL  = 0.8
#     p0  = 2e-8
#     d1  = 9
#     upcoeff = 12.2189
#     downcoeff = 0.9404
#     rho = 1e-3
#  [../]
#[]


[AuxKernels]
  [./accel_p]
    type = NewmarkAccelAux
    variable = accel_p
    displacement = p
    velocity = vel_p
    beta = 0.3025
    execute_on = timestep_end
  [../]
  [./vel_p]
    type = NewmarkVelAux
    variable = vel_p
    acceleration = accel_p
    gamma = 0.6
    execute_on = timestep_end
  [../]

[]


[BCs]
  [./top]
    type = CoupledNeumannBC
    variable = p
    boundary = top2 
    some_val = accel_z_in
    Reyold  = 1.0
  [../]

  [./inner]
    type = CoupledNeumannVectorBC
    variable = p
    boundary = inner
    some_val_x = accel_x_in
    some_val_y = accel_y_in
    Reyold  = -1.0
  [../]

  [./bottom]
    type = CoupledNeumannBC
    variable = p
    boundary = bottom2
    some_val = accel_z_in
    Reyold  = -1.0
  [../]
[]

[Materials]
  [./density]
    type = GenericConstantMaterial
    prop_names = 'density Diff'
    prop_values = '444.44 1000'
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
  l_tol = 1e-12
  nl_rel_tol = 1e-12
  dt = 5e-3
  l_max_its = 200

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
        file_base = AcousticOut
  [../]
[]

[MultiApps]
  [./sub]
    type = TransientMultiApp
    app_type = ASFracture
    execute_on = timestep_end
    positions = '0.0 0.0 1.0'
    input_files = asfractureDisp.i
  [../]
[]

[Transfers]
  [./from_sub_x]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = accel_x
    variable = accel_x_in
    source_boundary = 'inner'
    target_boundary = 'inner'
    execute_on = 'timestep_begin'
  [../]
  [./from_sub_y]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = accel_y
    variable = accel_y_in
    source_boundary = 'inner'
    target_boundary = 'inner'
    execute_on = 'timestep_begin'
  [../]
  [./from_sub_z_top]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = accel_z
    variable = accel_z_in
    source_boundary = 'top2'
    target_boundary = 'top2'
    execute_on = 'timestep_begin'
  [../]
  [./from_sub_z_bottom]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = accel_z
    variable = accel_z_in
    source_boundary = 'bottom2'
    target_boundary = 'bottom2'
    execute_on = 'timestep_begin'
  [../]

  [./to_sub_p_inner]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = p
    variable = pre_wave
    source_boundary = 'inner'
    target_boundary = 'inner'
    execute_on = 'timestep_end'
  [../]
  [./to_sub_p_top]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = p
    variable = pre_wave
    source_boundary = 'top2'
    target_boundary = 'top2'
    execute_on = 'timestep_end'
  [../]
  [./to_sub_p_bottom]
    type = MultiAppNearestNodeTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = p
    variable = pre_wave
    source_boundary = 'bottom2'
    target_boundary = 'bottom2'
    execute_on = 'timestep_end'
  [../]
[]
    
