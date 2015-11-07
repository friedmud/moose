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

#ifndef MULTICIRCLEEXACTGEOMETRY_H
#define MULTICIRCLEEXACTGEOMETRY_H

// Moose
#include "ExactGeometry.h"

class MultiCircleExactGeometry;

template<>
InputParameters validParams<MultiCircleExactGeometry>();

/**
 * Allows for popping new nodes generated in mesh adaptivity to circles.
 */
class MultiCircleExactGeometry :
  public ExactGeometry
{
public:
  MultiCircleExactGeometry(const InputParameters & parameters);
  virtual ~MultiCircleExactGeometry() {}

  /**
   * Move the passsed in Node to its correct position.
   *
   * Users should override this function and modify node to have the correct position.
   */
  virtual void moveNode(boundary_id_type boundary_id, Node & node);
};

#endif /* MULTICIRCLEEXACTGEOMETRY_H */
