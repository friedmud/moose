#include "RayTracingApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

// RayTracing Problems
#include "RayProblem.h"

// RayTracing Studies
#include "ElementCentroidRays.h"

// Ray Kernels
#include "AccumulateDistance.h"

// Ray Materials
#include "ConstantRayMaterial.h"

// Actions
#include "AddRayMaterialAction.h"
#include "AddRayKernelAction.h"
#include "AddRayBCAction.h"

template <>
InputParameters
validParams<RayTracingApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

RayTracingApp::RayTracingApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  RayTracingApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  RayTracingApp::associateSyntax(_syntax, _action_factory);
}

RayTracingApp::~RayTracingApp() {}

// External entry point for dynamic application loading
extern "C" void
RayTracingApp__registerApps()
{
  RayTracingApp::registerApps();
}
void
RayTracingApp::registerApps()
{
  registerApp(RayTracingApp);
}

// External entry point for dynamic object registration
extern "C" void
RayTracingApp__registerObjects(Factory & factory)
{
  RayTracingApp::registerObjects(factory);
}

void
RayTracingApp::registerObjects(Factory & factory)
{
  // RayTracing Problems
  registerProblem(RayProblem);

  // RayTracing Studies
  registerUserObject(ElementCentroidRays);

  // Ray Kernels
  registerUserObject(AccumulateDistance);

  // Ray Materials
  registerUserObject(ConstantRayMaterial);
}

// External entry point for dynamic syntax association
extern "C" void
RayTracingApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  RayTracingApp::associateSyntax(syntax, action_factory);
}
void
RayTracingApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  registerAction(AddRayKernelAction, "add_ray_kernel");
  syntax.registerActionSyntax("AddRayKernelAction", "RayKernels/*");
  registerMooseObjectTask("add_ray_kernel", RayKernel, false);

  addTaskDependency("add_ray_kernel", "add_kernel");

  registerAction(AddRayBoundaryConditionAction, "add_ray_boundary_condition");
  syntax.registerActionSyntax("AddRayBoundaryConditionAction", "RayBCs/*");
  registerMooseObjectTask("add_ray_boundary_condition", RayBoundaryCondition, false);

  addTaskDependency("add_ray_boundary_condition", "add_kernel");

  registerAction(AddRayMaterialAction, "add_ray_material");
  syntax.registerActionSyntax("AddRayMaterialAction", "RayMaterials/*");
  registerMooseObjectTask("add_ray_material", RayMaterial, false);

  addTaskDependency("add_ray_material", "add_kernel");
}
