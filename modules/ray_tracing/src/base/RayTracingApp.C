#include "RayTracingApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

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
RayTracingApp::registerObjects(Factory & /*factory*/)
{
}

// External entry point for dynamic syntax association
extern "C" void
RayTracingApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  RayTracingApp::associateSyntax(syntax, action_factory);
}
void
RayTracingApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
