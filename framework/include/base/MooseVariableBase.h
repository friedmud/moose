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

#ifndef MOOSEVARIABLEBASE_H
#define MOOSEVARIABLEBASE_H

#include "MooseTypes.h"
#include "MooseArray.h"

#include "libmesh/variable.h"
#include "libmesh/dof_map.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

// MetaPhysicL
#include "metaphysicl/dualnumber.h"
#include "metaphysicl/numberarray.h"

// The 100 here is for how many DoFs there are per element.
#define AD_MAX_DOFS_PER_ELEM 100
typedef MetaPhysicL::DualNumber<double, MetaPhysicL::NumberArray<AD_MAX_DOFS_PER_ELEM, double> > ADReal;

template<>
struct CompareTypes<double, ADReal>
{
  typedef ADReal supertype;
};

template<>
struct CompareTypes<ADReal, double>
{
  typedef ADReal supertype;
};

typedef VectorValue<ADReal> ADRealVectorValue;
typedef ADRealVectorValue ADRealGradient;

typedef TensorValue<ADReal> ADRealTensorValue;
typedef ADRealTensorValue ADRealTensor;

typedef MooseArray<Real>               VariableValue;
typedef MooseArray<RealGradient>       VariableGradient;
typedef MooseArray<RealTensor>         VariableSecond;

typedef MooseArray<ADReal>             ADVariableValue;
typedef MooseArray<ADRealGradient>     ADVariableGradient;
typedef MooseArray<ADRealTensor>       ADVariableSecond;

typedef MooseArray<std::vector<Real> >         VariableTestValue;
typedef MooseArray<std::vector<RealGradient> > VariableTestGradient;
typedef MooseArray<std::vector<RealTensor> >   VariableTestSecond;

typedef MooseArray<std::vector<Real> >         VariablePhiValue;
typedef MooseArray<std::vector<RealGradient> > VariablePhiGradient;
typedef MooseArray<std::vector<RealTensor> >   VariablePhiSecond;

class Assembly;
class SubProblem;
class SystemBase;


class MooseVariableBase
{
public:
  MooseVariableBase(unsigned int var_num, SystemBase & sys, Assembly & assembly, Moose::VarKindType var_kind);
  virtual ~MooseVariableBase();

  /**
   * Get variable number coming from libMesh
   * @return the libmesh variable number
   */
  unsigned int number() const { return _var_num; }

  /**
   * Get the system this variable is part of.
   */
  SystemBase & sys() { return _sys; }

  /**
   * Get the variable number
   */
  const std::string & name() const;

  /**
   * Kind of the variable (Nonlinear, Auxiliary, ...)
   */
  Moose::VarKindType kind() const { return _var_kind; }

  /**
   * Set the scaling factor for this variable
   */
  void scalingFactor(Real factor) { _scaling_factor = factor; }

  /**
   * Get the scaling factor for this variable
   */
  Real scalingFactor() const { return _scaling_factor; }

  /**
   * Get the order of this variable
   */
  unsigned int order() const;

  /**
   * The DofMap associated with the system this variable is in.
   */
  const DofMap & dofMap() { return _dof_map; }

  std::vector<dof_id_type> & dofIndices() { return _dof_indices; }

  unsigned int numberOfDofs() { return _dof_indices.size(); }

  /**
   * Is this variable nodal
   * @return true if it nodal, otherwise false
   */
  virtual bool isNodal() const = 0;

protected:
  /// variable number (from libMesh)
  unsigned int _var_num;
  /// variable number within MOOSE
  unsigned int _index;
  Moose::VarKindType _var_kind;
  /// Problem this variable is part of
  SubProblem & _subproblem;
  /// System this variable is part of
  SystemBase & _sys;

  /// libMesh variable object for this variable
  const Variable & _variable;

  /// Assembly data
  Assembly & _assembly;
  /// DOF map
  const DofMap & _dof_map;
  /// DOF indices
  std::vector<dof_id_type> _dof_indices;

  /// scaling factor for this variable
  Real _scaling_factor;
};

#endif /* MOOSEVARIABLEBASE_H */
