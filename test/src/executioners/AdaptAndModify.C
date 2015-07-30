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

#include "AdaptAndModify.h"
#include "TimeStepper.h"

//Moose includes

template<>
InputParameters validParams<AdaptAndModify>()
{
  InputParameters params = validParams<Transient>();
  params.addParam<unsigned int>("adapt_cycles", 1, "Number of adaptivity cycles to do.");
  return params;
}

AdaptAndModify::AdaptAndModify(const InputParameters & parameters) :
    Transient(parameters),
    _adapt_cycles(parameters.get<unsigned int>("adapt_cycles"))
{}

void
AdaptAndModify::incrementStepOrReject()
{
  if (lastSolveConverged())
  {
    _time_old = _time;
    _t_step++;

    _problem.advanceState();
  }
  else
  {
    _time_stepper->rejectStep();
    _time = _time_old;
  }

  _first = false;
}

void
AdaptAndModify::endTransientStep(Real input_time)
{
  Transient::endTransientStep(input_time);

  if (lastSolveConverged())
  {
    // Compute the Error Indicators and Markers
    for (unsigned int i=0; i<_adapt_cycles; i++)
    {
      // Compute the Error Indicators and Markers
      _problem.computeIndicatorsAndMarkers();

#ifdef LIBMESH_ENABLE_AMR
      if (_problem.adaptivity().isOn())
        _problem.adaptMesh();

#endif
    }
  }

}
