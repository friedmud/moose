#ifndef RAYTRACINGAPP_H
#define RAYTRACINGAPP_H

#include "MooseApp.h"

class RayTracingApp;

template <>
InputParameters validParams<RayTracingApp>();

class RayTracingApp : public MooseApp
{
public:
  RayTracingApp(InputParameters parameters);
  virtual ~RayTracingApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* RAYTRACINGAPP_H */
