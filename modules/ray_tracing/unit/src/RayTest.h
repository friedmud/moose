/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef RAYTEST_H
#define RAYTEST_H

#include "gtest/gtest.h"

// Moose includes
#include "InputParameters.h"
#include "FEProblem.h"
#include "RayTracingApp.h"
#include "AppFactory.h"
#include "GeneratedMesh.h"

// libMesh includes
#include "libmesh/plane.h"

// System Includes
#include <memory>

class RayTest : public ::testing::Test
{
protected:
  void SetUp()
  {
    const char * argv[2] = {"foo", "\0"};

    _app.reset(AppFactory::createApp("RayTracingApp", 1, (char **)argv));
    _factory = &_app->getFactory();

    InputParameters mesh_params = _factory->getValidParams("GeneratedMesh");
    mesh_params.set<MooseEnum>("dim") = "2";
    mesh_params.set<std::string>("_object_name") = "mesh";
    _mesh = libmesh_make_unique<GeneratedMesh>(mesh_params);
    _mesh->buildMesh();
  }

  std::unique_ptr<MooseApp> _app;
  Factory * _factory;
  std::unique_ptr<MooseMesh> _mesh;
};

#endif // RAYTEST_H
