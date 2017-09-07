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

#ifndef RAYKERNEL_H
#define RAYKERNEL_H

// Local Includes
#include "Ray.h"

// MOOSE Includes
#include "MooseObject.h"
#include "SetupInterface.h"
#include "Restartable.h"
#include "BlockRestrictable.h"
#include "Coupleable.h"

class RayKernel;
class RayProblemBase;
class RaySystem;

template <>
InputParameters validParams<RayKernel>();

class RayKernel : public MooseObject,
                  public SetupInterface,
                  public Restartable,
                  public BlockRestrictable,
                  public Coupleable
{
public:
  RayKernel(const InputParameters & params);
  virtual ~RayKernel();

  virtual void initialSetup() {}
  virtual void timestepSetup() {}

  /**
   * Called at the beginning of a new Ray trace
   *
   * Useful for caching data!
   */
  virtual void rayStart() {}

  /**
   * Called on each Segment
   * @param start The beginning of the segment
   * @param end The end of the segment
   * @param ends_in_elem Whether or not the Ray ends _within_ the current element.  Note: this is
   * NOT true if the Ray hits a physical boundary on one side of the element.  It is _only_ true if
   * the Ray truly ends _within_ the element.
   */
  virtual void onSegment(const Elem * /*elem*/,
                         const Point & /*start*/,
                         const Point & /*end*/,
                         bool /*ends_in_elem*/) = 0;

  /**
   * Set the current Ray that's being worked on
   */
  virtual void setRay(const std::shared_ptr<Ray> & ray) { _ray = ray.get(); }

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

  /// Number of groups
  unsigned int _num_groups;

  /// The Ray that's being worked on
  Ray * _ray;

  /// Number of polar angles
  unsigned int _num_polar;

  /// Offest into the vectors associated with the RaySystem
  dof_id_type & _current_offset;

  /// The current group values for this thread
  PetscScalar *& _group_solution_values;
};

#endif /* RAYKERNEL_H */
