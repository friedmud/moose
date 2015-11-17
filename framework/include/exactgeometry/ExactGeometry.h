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

#ifndef EXACTGEOMETRY_H
#define EXACTGEOMETRY_H

// Moose
#include "MooseObject.h"
#include "InputParameters.h"
#include "SetupInterface.h"
#include "BoundaryRestrictable.h"
#include "Restartable.h"

class ExactGeometry;
class SubProblem;
class FEProblem;

template<>
InputParameters validParams<ExactGeometry>();

/**
 * Base class for all ExactGeometry objects.
 *
 * ExactGeometry allows for moving new nodes created in
 * mesh adaptivity to the real geometry.
 */
class ExactGeometry :
  public MooseObject,
  public SetupInterface,
  public BoundaryRestrictable,
  public Restartable
{
public:
  ExactGeometry(const InputParameters & parameters);
  virtual ~ExactGeometry() {}

  /**
   * Move the passsed in Node to its correct position.
   *
   * Users should override this function and modify node to have the correct position.
   */
  virtual void moveNode(boundary_id_type boundary_id, Node & node) = 0;

protected:
  SubProblem & _subproblem;
  FEProblem & _fe_problem;

  THREAD_ID _tid;
};

#endif /* EXACTGEOMETRY_H */
