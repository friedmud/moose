//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RegisterConstraints.h"

#include "Factory.h"

#include "TiedValueConstraint.h"
#include "CoupledTiedValueConstraint.h"
#include "AddBoundsVectorsAction.h"
#include "EqualGradientConstraint.h"
#include "EqualValueConstraint.h"
#include "EqualValueBoundaryConstraint.h"
#include "LinearNodalConstraint.h"

namespace Moose
{

/**
 * Called from Moose.C to register Actions
 */
void
registerConstraints(Factory & factory)
{
  registerConstraint(TiedValueConstraint);
  registerConstraint(CoupledTiedValueConstraint);
  registerConstraint(EqualGradientConstraint);
  registerConstraint(EqualValueConstraint);
  registerConstraint(EqualValueBoundaryConstraint);
  registerConstraint(LinearNodalConstraint);
}
}
