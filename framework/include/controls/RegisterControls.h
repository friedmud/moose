//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef REGISTER_CONTROLS_H
#define REGISTER_CONTROLS_H

// Forward Declarations
class Factory;

namespace Moose
{
/**
 * Called from Moose.C to register Controls
 */
void registerControls(Factory & factory);
}

#endif
