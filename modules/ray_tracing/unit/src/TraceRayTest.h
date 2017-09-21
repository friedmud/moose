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

#ifndef TRACERAYTEST_H
#define TRACERAYTEST_H

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

class TraceRayTest : public ::testing::Test
{
protected:
  void SetUp()
  {
    const char * argv[2] = {"foo", "\0"};

    _app.reset(AppFactory::createApp("RayTracingApp", 1, (char **)argv));
    _factory = &_app->getFactory();
  }

  std::unique_ptr<MooseApp> _app;
  Factory * _factory;
};

#endif // TRACERAYTEST_H
