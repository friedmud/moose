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

#ifndef CONSTANTGROUPVALUES_H
#define CONSTANTGROUPVALUES_H

// Local Includes
#include "Ray.h"
#include "RayKernel.h"

class ConstantGroupValues;

template <>
InputParameters validParams<ConstantGroupValues>();

class ConstantGroupValues : public RayKernel
{
public:
  ConstantGroupValues(const InputParameters & params);
  virtual ~ConstantGroupValues();

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
                         bool /*ends_in_elem*/);

  /**
   * Set the current Ray that's being worked on
   */
  virtual void setRay(const std::shared_ptr<Ray> & ray);

protected:
  /// Number of groups in the problem
  unsigned int _num_groups;

  /// Values to use for each group
  const std::vector<Real> & _input_group_values;

  /// Pointer to the beginning of the Ray's data
  Real * _ray_data;
};

#endif /* CONSTANTGROUPVALUES_H */
