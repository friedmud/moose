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

#ifndef AVERAGEGROUPVALUEAUX_H
#define AVERAGEGROUPVALUEAUX_H

#include "AuxKernel.h"

// Local Includes
#include "GroupValueAuxBase.h"

// Forward Declarations
class AverageGroupValueAux;

template <>
InputParameters validParams<AverageGroupValueAux>();

class AverageGroupValueAux : public GroupValueAuxBase
{
public:
  AverageGroupValueAux(const InputParameters & parameters);

  virtual ~AverageGroupValueAux() {}

protected:
  virtual Real computeValue();
};

#endif // AVERAGEGROUPVALUEAUX_H
