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

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
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
  num_groups = 1
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
  [./ad]
    type = AccumulateDistance
  [../]
[]

[RayMaterials]
  [./crm]
    type = ConstantRayMaterial
    sigma_t = 1
  [../]
[]

