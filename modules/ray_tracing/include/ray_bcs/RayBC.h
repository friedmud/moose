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

#ifndef RAYBC_H
#define RAYBC_H

// Local Includes
#include "Ray.h"

// MOOSE Includes
#include "MooseObject.h"
#include "SetupInterface.h"
#include "Restartable.h"
#include "BoundaryRestrictableRequired.h"
#include "Coupleable.h"
#include "GeometricSearchInterface.h"

class RayBC;
class RayProblemBase;
class RaySystem;

template <>
InputParameters validParams<RayBC>();

class RayBC : public MooseObject,
              public SetupInterface,
              public Restartable,
              public BoundaryRestrictableRequired,
              public Coupleable,
              public GeometricSearchInterface
{
public:
  RayBC(const InputParameters & params);
  virtual ~RayBC();

  virtual void initialSetup() {}
  virtual void timestepSetup() {}

  /**
   * Called on the boundary.
   *
   * Must modify the Ray
   */
  virtual void apply(const Elem * elem,
                     const unsigned int intersected_side,
                     const Point & intersection_point,
                     const std::shared_ptr<Ray> & ray) = 0;

  /**
   * These shouldn't really be used
   */
  virtual void initialize(){};
  virtual void execute(){};
  virtual void finalize(){};

  /**
   * This needs to be here to satisfy the MooseObjectWarehouse interface
   */
  class DummyVariable
  {
  public:
    std::string name() { return ""; }
  };

  /**
   * This needs to be here to satisfy the MooseObjectWarehouse interface
   */
  DummyVariable variable() { return DummyVariable(); }

protected:
  /// The Ray Problem
  RayProblemBase & _ray_problem;

  /// The Ray System
  RaySystem & _ray_sys;

  /// The thread ID this object is assigned to
  THREAD_ID _tid;
};

#endif /* RAYBC_H */
