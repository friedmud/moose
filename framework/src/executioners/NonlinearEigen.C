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

#include "NonlinearEigen.h"
#include "EigenSystem.h"
#include "Function.h"

template<>
InputParameters validParams<NonlinearEigen>()
{
  InputParameters params = validParams<EigenExecutionerBase>();
  params.addParam<unsigned int>("free_power_iterations", 4, "The number of free power iterations");
  params.addParam<Real>("source_abs_tol", 1e-06, "Absolute tolernance on residual norm");
  params.addParam<Real>("source_rel_tol", 1e-50, "Relative tolernance on residual norm after free power iterations");
  params.addParam<Real>("pfactor", 1e-2, "The factor of residual to be reduced per power iteration");
  params.addParam<Real>("k0", 1.0, "Initial guess of the eigenvalue");
  params.addParam<bool>("output_pi_history", false, "True to output solutions durint PI");
  params.addRequiredParam<FunctionName>("power_modulating_function", "A Function of time with values between 0 and 1 that will 'modulate' the max power over time");
  return params;
}

NonlinearEigen::NonlinearEigen(const std::string & name, InputParameters parameters)
    :EigenExecutionerBase(name, parameters),
     // local static memebers
     _free_iter(getParam<unsigned int>("free_power_iterations")),
     _abs_tol(getParam<Real>("source_abs_tol")),
     _rel_tol(getParam<Real>("source_rel_tol")),
     _pfactor(getParam<Real>("pfactor")),
     _output_pi(getParam<bool>("output_pi_history")),
     _power_modulating_function(NULL)
{
  _eigenvalue = getParam<Real>("k0");
  addRealParameterReporter("eigenvalue");

  if (getParam<bool>("output_on_final") && _output_pi)
  {
    mooseWarning("Only final solution will be outputted, output_pi_history=true will be ignored!");
    _output_pi = false;
  }
}

void
NonlinearEigen::init()
{
  EigenExecutionerBase::init();

  if (_free_iter>0)
  {
    // save the initial guess
    _problem.advanceState();

    // free power iterations
    Moose::out << std::endl << " Free power iteration starts"  << std::endl;

    Real initial_res;
    inversePowerIteration(_free_iter, _free_iter, _pfactor, false,
                          std::numeric_limits<Real>::min(), std::numeric_limits<Real>::max(),
                          true, _output_pi, 0.0, _eigenvalue, initial_res);

    _problem.computeUserObjects(EXEC_TIMESTEP, UserObjectWarehouse::PRE_AUX);
    _problem.onTimestepEnd();
    _problem.computeAuxiliaryKernels(EXEC_TIMESTEP);
    _problem.computeUserObjects(EXEC_TIMESTEP, UserObjectWarehouse::POST_AUX);

    if (!getParam<bool>("output_on_final"))
    {
      _problem.timeStep()++;
      Real t = _problem.time();
      _problem.time() = _problem.timeStep();
      _output_warehouse.outputStep();
      _problem.time() = t;
    }
  }
}

void
NonlinearEigen::execute()
{
  _power_modulating_function = &_problem.getFunction(getParam<FunctionName>("power_modulating_function"));

  preExecute();

  for(unsigned int step=0; step < 40; step++)
  {
    Real multi_app_dt = _problem.computeMultiAppsDT(EXEC_TIMESTEP_BEGIN);
    _problem.dt() = multi_app_dt; // This means that we are going

    multi_app_dt = _problem.computeMultiAppsDT(EXEC_TIMESTEP);
    if(multi_app_dt < _problem.dt())
      _problem.dt() = multi_app_dt;

    //_problem.dt() *= _dt_modulating_function->value(_problem.time(), Point());

    _problem.time() += _problem.dt();

    takeStep();

    Real scaling = _power_modulating_function->value(_problem.time(), Point()) * getParam<Real>("normal_factor");

    // Scale the solution to get the correct power
    _eigen_sys.scaleSystemSolution(EigenSystem::EIGEN, scaling);
    Moose::out << " Solution is rescaled with factor " << scaling << " for normalization!" << std::endl;


    // update all aux variables and user objects
    for (unsigned int i=0; i<Moose::exec_types.size(); i++)
    {
      // EXEC_CUSTOM is special, should be treated only by specifically designed executioners.
      if (Moose::exec_types[i]==EXEC_CUSTOM) continue;
      _problem.computeUserObjects(Moose::exec_types[i], UserObjectWarehouse::PRE_AUX);
      _problem.computeAuxiliaryKernels(Moose::exec_types[i]);
      _problem.computeUserObjects(Moose::exec_types[i], UserObjectWarehouse::POST_AUX);
    }


    _output_warehouse.outputStep();

    // Solve BISON
    _problem.execTransfers(EXEC_TIMESTEP_BEGIN);
    _problem.execMultiApps(EXEC_TIMESTEP_BEGIN);

    // Solve R7
    _problem.execTransfers(EXEC_TIMESTEP);
    _problem.execMultiApps(EXEC_TIMESTEP);

    // Scale back
    _eigen_sys.scaleSystemSolution(EigenSystem::EIGEN, 1.0/scaling);

//    postExecute();
  }
}

void
NonlinearEigen::takeStep()
{
  Moose::out << " Nonlinear iteration starts"  << std::endl;

  // nonlinear solve
  _problem.computeUserObjects(EXEC_TIMESTEP_BEGIN, UserObjectWarehouse::PRE_AUX);
  preSolve();
  _problem.timestepSetup();
  _problem.advanceState();
  _problem.computeUserObjects(EXEC_TIMESTEP_BEGIN, UserObjectWarehouse::POST_AUX);

  nonlinearSolve(_rel_tol, _abs_tol, _pfactor, _eigenvalue);

  postSolve();
  printEigenvalue();

  _problem.computeUserObjects(EXEC_TIMESTEP, UserObjectWarehouse::PRE_AUX);
  _problem.onTimestepEnd();
  _problem.computeAuxiliaryKernels(EXEC_TIMESTEP);
  _problem.computeUserObjects(EXEC_TIMESTEP, UserObjectWarehouse::POST_AUX);
}
