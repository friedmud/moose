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

#include "ConstantRayMaterial.h"

// Local Includes
#include "RayProblem.h"

// MOOSE Includes
#include "MooseMesh.h"

template <>
InputParameters
validParams<ConstantRayMaterial>()
{
  InputParameters params = validParams<RayMaterial>();

  params.addParam<std::vector<Real>>("sigma_t", "Total");

  return params;
}

ConstantRayMaterial::ConstantRayMaterial(const InputParameters & params)
  : RayMaterial(params), _sigma_t(getParam<std::vector<Real>>("sigma_t"))
{
}
