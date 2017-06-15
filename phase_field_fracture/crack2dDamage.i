#Fix y displacement on the left side
[Mesh]
  type = FileMesh
  file = crack_mesh.e
  uniform_refine = 1
[]

#[GlobalParams]
#  displacements = 'disp_x disp_y'
#[]

[Variables]
  [./d]
  [../]
  [./b]
  [../]
[]

[AuxVariables]
  [./disp_x_in]
  [../]
  [./disp_y_in]
  [../]
[]

[Kernels]
  [./pfbulk]
     type = CohesivePFFracBulkRate
     variable = d
     ifOld = true
     l = 0.05
     p = 3
     beta = b
     visco = 1.e-4
     gc_prop_var = 'gc_prop'
     G0_var = 'G0_pos'
     dG0_dstrain_var = 'dG0_pos_dstrain'
     Emod = 'E'
     sigmac = 'sc'
     disp_x = disp_x_in
     disp_y = disp_y_in
  [../]
  [./dcdt]
     type = TimeDerivativeExp
     variable = d
#      lumping = true
  [../]
  [./pfintvar]
      type = PFFracIntVar
      variable = b
  [../]
  [./pfintcoupled]
      type = PFFracCoupledInterfaceExp
      variable = b
      c = d
   [../]
[]


[Materials]
  [./pfbulkmat]
    type = PFFracBulkRateMaterial
    gc = 2.7e-3
  [../]

  [./elastic]
       type = CohesiveLinearIsoElasticPFDamage
       c = d_in
       kdamage = 1e-6
       store_stress_old = true
       gc_prop_var = 'gc_prop'
       Emod = 'E'
       sigmac = 'sc'
       l = 0.05
       p = 3
  [../]

  [./constant]
    type = GenericConstantMaterial
    prop_names = 'E sc'
    prop_values = '210.0 0.8646'
  [../]

  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '121.0 81.0'
    fill_method = symmetric_isotropic
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x_in disp_y_in'
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
#  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
#  petsc_options_value = 'asm      31                  preonly       lu           1'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  nl_rel_tol = 1e-8
  l_max_its = 30
  nl_max_its = 30

  dt = 1e-4
  dtmin = 1e-10
  start_time = 0.0
  end_time = 1
  #num_steps = 2
[]

[Outputs]
  file_base = ModeIIDamage
  interval = 10
  exodus = true
  csv = true
  gnuplot = true
[]
