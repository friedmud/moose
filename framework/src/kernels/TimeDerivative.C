//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TimeDerivative.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableField.h"

#include "libmesh/quadrature.h"

registerMooseObject("MooseApp", TimeDerivative);

template <>
InputParameters
validParams<TimeDerivative>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addClassDescription("The time derivative operator with the weak form of $(\\psi_i, "
                             "\\frac{\\partial u_h}{\\partial t})$.");
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");

  return params;
}

TimeDerivative::TimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters), _lumping(getParam<bool>("lumping"))
{
}

Real
TimeDerivative::computeQpResidual()
{
  return _test[_i][_qp] * _u_dot[_qp];
}

Real
TimeDerivative::computeQpJacobian()
{
  return _test[_i][_qp] * _phi[_j][_qp] * _du_dot_du[_qp];
}

void
TimeDerivative::computeJacobian()
{
  if (_lumping)
  {
    prepareMatrixTag(_assembly, _var.number(), _var.number());

    precalculateJacobian();
    auto test_size = _test.size();
    auto phi_size = _phi.size();
    auto n_points = _qrule->n_points();

    for (_i = 0; _i < test_size; _i++)
      for (_j = 0; _j < phi_size; _j++)
        for (_qp = 0; _qp < n_points; _qp++)
          _local_ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();

    accumulateTaggedLocalMatrix();
  }
  else
    TimeKernel::computeJacobian();
}
