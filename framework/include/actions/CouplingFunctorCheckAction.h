//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"

class CouplingFunctorCheckAction;

template <>
InputParameters validParams<CouplingFunctorCheckAction>();

/**
 * Checks whether there are any Kernels or BoundaryConditions in the warehouses and if so adds a
 * default coupling functor to ensure correct sparsity
 */
class CouplingFunctorCheckAction : public Action
{
public:
  static InputParameters validParams();

  CouplingFunctorCheckAction(InputParameters parameters);

protected:
  void act() override;

  const PerfID _reinitializing_vectors_because_of_algebraic_ghosting_timer;
  const PerfID _adding_relationship_managers_timer;
  const PerfID _reiniting_nonlinear_system;
};
