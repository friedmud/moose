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
  [./rays_started]
    family = MONOMIAL
    order = CONSTANT
  []

  [./rays_traced]
    family = MONOMIAL
    order = CONSTANT
  []

  [./chunks_traced]
    family = MONOMIAL
    order = CONSTANT
  []

  [./rays_received]
    family = MONOMIAL
    order = CONSTANT
  []

  [./buffers_received]
    family = MONOMIAL
    order = CONSTANT
  []

  [./rays_sent]
    family = MONOMIAL
    order = CONSTANT
  []

  [./buffers_sent]
    family = MONOMIAL
    order = CONSTANT
  []

  [./intersections]
    family = MONOMIAL
    order = CONSTANT
  []

  [./generation_time]
    family = MONOMIAL
    order = CONSTANT
  []

  [./propagation_time]
    family = MONOMIAL
    order = CONSTANT
  []

  [./ray_pool_created]
    family = MONOMIAL
    order = CONSTANT
  []

  [./receive_ray_pool_created]
    family = MONOMIAL
    order = CONSTANT
  []

  [./receive_buffer_pool_created]
    family = MONOMIAL
    order = CONSTANT
  []

  [./send_buffer_pool_created]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[AuxKernels]
  #These are inactive because they aren't stable to gold against
  inactive = 'generation_time propagation_time'

  [./rays_started]
    type = VectorPostprocessorVisualizationAux
    variable = rays_started
    vpp = per_proc_ray_tracing
    vector_name = rays_started
    execute_on = timestep_end
  []

  [./rays_traced]
    type = VectorPostprocessorVisualizationAux
    variable = rays_traced
    vpp = per_proc_ray_tracing
    vector_name = rays_traced
    execute_on = timestep_end
  []

  [./chunks_traced]
    type = VectorPostprocessorVisualizationAux
    variable = chunks_traced
    vpp = per_proc_ray_tracing
    vector_name = chunks_traced
    execute_on = timestep_end
  []

  [./rays_received]
    type = VectorPostprocessorVisualizationAux
    variable = rays_received
    vpp = per_proc_ray_tracing
    vector_name = rays_received
    execute_on = timestep_end
  []

  [./buffers_received]
    type = VectorPostprocessorVisualizationAux
    variable = buffers_received
    vpp = per_proc_ray_tracing
    vector_name = buffers_received
    execute_on = timestep_end
  []

  [./rays_sent]
    type = VectorPostprocessorVisualizationAux
    variable = rays_sent
    vpp = per_proc_ray_tracing
    vector_name = rays_sent
    execute_on = timestep_end
  []

  [./buffers_sent]
    type = VectorPostprocessorVisualizationAux
    variable = buffers_sent
    vpp = per_proc_ray_tracing
    vector_name = buffers_sent
    execute_on = timestep_end
  []

  [./intersections]
    type = VectorPostprocessorVisualizationAux
    variable = intersections
    vpp = per_proc_ray_tracing
    vector_name = intersections
    execute_on = timestep_end
  []

  [./generation_time]
    type = VectorPostprocessorVisualizationAux
    variable = generation_time
    vpp = per_proc_ray_tracing
    vector_name = generation_time
    execute_on = timestep_end
  []

  [./propagation_time]
    type = VectorPostprocessorVisualizationAux
    variable = propagation_time
    vpp = per_proc_ray_tracing
    vector_name = propagation_time
    execute_on = timestep_end
  []

  [./ray_pool_created]
    type = VectorPostprocessorVisualizationAux
    variable = ray_pool_created
    vpp = per_proc_ray_tracing
    vector_name = ray_pool_created
    execute_on = timestep_end
  []

  [./receive_ray_pool_created]
    type = VectorPostprocessorVisualizationAux
    variable = receive_ray_pool_created
    vpp = per_proc_ray_tracing
    vector_name = receive_ray_pool_created
    execute_on = timestep_end
  []

  [./receive_buffer_pool_created]
    type = VectorPostprocessorVisualizationAux
    variable = receive_buffer_pool_created
    vpp = per_proc_ray_tracing
    vector_name = receive_buffer_pool_created
    execute_on = timestep_end
  []

  [./send_buffer_pool_created]
    type = VectorPostprocessorVisualizationAux
    variable = send_buffer_pool_created
    vpp = per_proc_ray_tracing
    vector_name = send_buffer_pool_created
    execute_on = timestep_end
  []
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
    ray_distance = 10000
    num_rays = 10000
    end_elems = '3 4'
    start_elems = '0 1'
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

[VectorPostprocessors]
  [./per_proc_ray_tracing]
    type = PerProcessorRayTracingResultsVectorPostprocessor
    results = 'rays_started rays_traced chunks_traced rays_received buffers_received rays_sent buffers_sent intersections generation_time propagation_time ray_pool_created receive_ray_pool_created receive_buffer_pool_created send_buffer_pool_created'
    execute_on = 'timestep_end'
  [../]
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

[Outputs]
  exodus = true
[]
