[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 5
  xmax = 5
[]

[Variables]
  [./u]
  [../]
[]

[AuxVariables]
  [./group_average]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[AuxKernels]
  [./average_group_value]
    type = AverageGroupValueAux
    variable = group_average
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]
  [./right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
[]

[UserObjects]
  [./ecr]
    type = ElementCentroidRays
    ray_distance = 1e5
    num_rays = 1e5
    end_elems = 4
    start_elems = 0
  [../]
[]

[Problem]
  type = RayProblem
  num_groups = 3
  study = ecr
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]

[RayKernels]
  [./cgv]
    group_values = '1 2 3'
    type = ConstantGroupValues
  [../]
[]

[RayMaterials]
  [./crm]
    type = ConstantRayMaterial
    sigma_t = 1
  [../]
[]

