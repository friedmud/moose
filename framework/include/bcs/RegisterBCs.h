//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef REGISTER_BCS_H
#define REGISTER_BCS_H

// Forward Declarations
class Syntax;
class ActionFactory;

namespace Moose
{
/**
 * Called from Moose.C to register BCs
 */
void registerBCs(Factory & factory);
}

#endif
