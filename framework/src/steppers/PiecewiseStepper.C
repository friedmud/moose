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

#include "PiecewiseStepper.h"
#include "FEProblem.h"
#include "Transient.h"
#include "MooseApp.h"

template<>
InputParameters validParams<PiecewiseStepper>()
{
  InputParameters params = validParams<Stepper>();

  params.addRequiredParam<std::vector<Real> >("times", "The values of t");
  params.addRequiredParam<std::vector<Real> >("dts",   "The values of dt");
  params.addParam<bool>("interpolate", true, "Whether or not to interpolate DT between times.  This is true by default for historical reasons.");

  return params;
}

PiecewiseStepper::PiecewiseStepper(const InputParameters & parameters) :
    Stepper(parameters),
    _times(getParam<std::vector<Real> >("times")),
    _dts(getParam<std::vector<Real> >("dts")),
    _interpolate(getParam<bool>("interpolate")),
    _linear_interpolation(_times, _dts)
{
}

Real
PiecewiseStepper::computeDT()
{
  if (_interpolate)
    return _linear_interpolation.sample(_time);

  if (MooseUtils::relativeFuzzyGreaterEqual(_time, _times.back()))
    return _dts.back();

  for (int i = 0; i < _times.size() - 1; i++)
    if (MooseUtils::relativeFuzzyLessThan(_time, _times[i + 1]))
      return _dts[i];

  return _dts.back();
}

Real
PiecewiseStepper::computeFailedDT()
{
  return 0.5 * _dt[0]; // _dt[0] was the dt actually used in the last timestep
}

PiecewiseStepper::~PiecewiseStepper()
{
}
