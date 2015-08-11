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

#include "Steady.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseApp.h"
#include "libmesh/equation_systems.h"

template<>
InputParameters validParams<Steady>()
{
  return validParams<Executioner>();
}


Steady::Steady(const InputParameters & parameters) :
    Executioner(parameters)
{
  _fe_problem.getNonlinearSystem().setDecomposition(_splitting);

  if (!_restart_file_base.empty())
    _fe_problem.setRestartFile(_restart_file_base);
  {
    std::string ti_str = "SteadyState";
    InputParameters params = _app.getFactory().getValidParams(ti_str);
    _fe_problem.addTimeIntegrator(ti_str, "ti", params);
  }
}

Steady::~Steady()
{
}

void
Steady::_init()
{
  if (_app.isRecovering())
  {
    _console << "\nCannot recover steady solves!\nExiting...\n" << std::endl;
    return;
  }

  // Everything else defaults to 1
  setCycles(_fe_problem.adaptivity().getSteps()+1);

  std::cout<<"Cycles: "<<_fe_problem.adaptivity().getSteps()+1<<std::endl;


  setTimeScheme(PSEUDO_TIME);
}

bool
Steady::keepStepping()
{
  if (_app.isRecovering())
    return false;

  return true;
}

bool
Steady::keepCycling()
{
  if (currentCycle() > 1)
    return _fe_problem.converged();
  else
    return true;
}

void
Steady::beginStep()
{
}

void
Steady::execute()
{
  _fe_problem.solve();
}

void
Steady::checkIntegrity()
{
  // check to make sure that we don't have any time kernels in this simulation (Steady State)
  if (_fe_problem.getNonlinearSystem().containsTimeKernel())
    mooseError("You have specified time kernels in your steady state simulation");
}
