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

#ifndef GROUPVALUEAUXBASE_H
#define GROUPVALUEAUXBASE_H

#include "AuxKernel.h"

// Local Includes
#include "RayProblem.h"

// Forward Declarations
class GroupValueAuxBase;

template <>
InputParameters validParams<GroupValueAuxBase>();

class GroupValueAuxBase : public AuxKernel
{
public:
  GroupValueAuxBase(const InputParameters & parameters);

  virtual ~GroupValueAuxBase() {}

protected:
  /// The RAY Problem
  RayProblem & _ray_problem;

  /// The RAY System
  RaySystem & _ray_sys;

  /// Number of energy groups in the problem
  unsigned int _num_groups;

  /// Offest into the vectors
  dof_id_type & _current_offset;

  /// READ only! The current group values
  PetscScalar *& _group_values;
};

#endif // GROUPVALUEAUXBASE_H
