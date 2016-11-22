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

#ifndef NEWSTEPPERINTERFACE_H
#define NEWSTEPPERINTERFACE_H

// Standard includes
#include <string>

// MOOSE includes
#include "MooseTypes.h"

// Forward Declarations
class FEProblem;
class InputParameters;
class StepperName;
class MooseObject;

/**
 * Interface class for classes which interact with Postprocessors.
 * Provides the getPostprocessorValueXYZ() and related interfaces.
 */
class NewStepperInterface
{
public:
  NewStepperInterface(const MooseObject * moose_object);

  /**
   * Get the DT from another Stepper.
   *
   * Note: This is intended to be called in the constructor of a Stepper
   * This returns a _reference_ and must be caught as a reference.
   * The value within this reference will automatically be updated by the system.
   *
   * @param name The name of the InputParameter holding the Stepper name to get the value from.
   */
  const Real & getStepperDT(const std::string & name);


  /**
   * Get the DT from another Stepper.
   *
   * Note: This is intended to be called in the constructor of a Stepper
   * This returns a _reference_ and must be caught as a reference.
   * The value within this reference will automatically be updated by the system.
   *
   * @param name The name of the Stepper to get the value from.
   */
  const Real & getStepperDTByName(const StepperName & name);

protected:
  /**
   * Check to see if a Stepper exists
   *
   * @param name The name of the parameter holding the name of the Stepper
   */
  bool hasStepper(const std::string & name) const;

  /**
   * Check to see if a Stepper exists
   *
   * @param name The name of the Stepper
   */
  bool hasStepperByName(const PostprocessorName & name) const;

private:
  /// NewStepperInterface Parameters
  const InputParameters & _si_params;

  /// Reference the the FEProblem class
  FEProblem & _si_feproblem;
};

#endif //NEWSTEPPERINTERFACE_H
