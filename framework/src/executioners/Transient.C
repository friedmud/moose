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

#include "Transient.h"

//Moose includes
#include "Factory.h"
#include "SubProblem.h"
#include "TimePeriod.h"
#include "TimeStepper.h"
#include "MooseApp.h"
#include "Conversion.h"
//libMesh includes
#include "libmesh/implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/numeric_vector.h"

// C++ Includes
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

template<>
InputParameters validParams<Transient>()
{
  InputParameters params = validParams<Executioner>();
  std::vector<Real> sync_times(1);
  sync_times[0] = -std::numeric_limits<Real>::max();

  /**
   * For backwards compatibility we'll allow users to set the TimeIntegration scheme inside of the executioner block
   * as long as the TimeIntegrator does not have any additional parameters.
   */
  MooseEnum schemes("implicit-euler explicit-euler crank-nicolson bdf2 rk-2 dirk");

  params.addParam<Real>("start_time",      0.0,    "The start time of the simulation");
  params.addParam<Real>("end_time",        1.0e30, "The end time of the simulation");
  params.addParam<Real>("dt",              1.,     "The timestep size between solves");
  params.addParam<Real>("dtmin",           2.0e-14,    "The minimum timestep size in an adaptive run");
  params.addParam<Real>("dtmax",           1.0e30, "The maximum timestep size in an adaptive run");
  params.addParam<bool>("reset_dt", false, "Use when restarting a calculation to force a change in dt.");
  params.addParam<unsigned int>("num_steps",       std::numeric_limits<unsigned int>::max(),     "The number of timesteps in a transient run");
  params.addParam<int> ("n_startup_steps", 0,      "The number of timesteps during startup");
  params.addParam<bool>("trans_ss_check",  false,  "Whether or not to check for steady state conditions");
  params.addParam<Real>("ss_check_tol",    1.0e-08,"Whenever the relative residual changes by less than this the solution will be considered to be at steady state.");
  params.addParam<Real>("ss_tmin",         0.0,    "Minimum number of timesteps to take before checking for steady state conditions.");
  params.addParam<Real>("predictor_scale", "The scale factor for the predictor (can range from 0 to 1)");

  params.addParam<std::vector<std::string> >("time_periods", "The names of periods");
  params.addParam<std::vector<Real> >("time_period_starts", "The start times of time periods");
  params.addParam<std::vector<Real> >("time_period_ends", "The end times of time periods");
  params.addParam<bool>("abort_on_solve_fail", false, "abort if solve not converged rather than cut timestep");
  params.addParam<MooseEnum>("scheme",          schemes,  "Time integration scheme used.");
  params.addParam<Real>("timestep_tolerance", 2.0e-14, "the tolerance setting for final timestep size and sync times");

  params.addParam<bool>("use_multiapp_dt", false, "If true then the dt for the simulation will be chosen by the MultiApps.  If false (the default) then the minimum over the master dt and the MultiApps is used");

  params.addParamNamesToGroup("start_time dtmin dtmax n_startup_steps trans_ss_check ss_check_tol ss_tmin sync_times time_t time_dt growth_factor predictor_scale use_AB2 use_littlef abort_on_solve_fail output_to_file file_name estimate_time_error timestep_tolerance use_multiapp_dt", "Advanced");

  params.addParamNamesToGroup("time_periods time_period_starts time_period_ends", "Time Periods");


  params.addParam<bool>("verbose", false, "Print detailed diagnostics on timestep calculation");

  return params;
}

Transient::Transient(const InputParameters & parameters) :
    Executioner(parameters),
    _problem(*parameters.getCheckedPointerParam<FEProblem *>("_fe_problem", "This might happen if you don't have a mesh")),
    _time_scheme(getParam<MooseEnum>("scheme")),
    _t_step(_problem.timeStep()),
    _time(_problem.time()),
    _time_old(_problem.timeOld()),
    _dt(_problem.dt()),
    _dt_old(_problem.dtOld()),
    _unconstrained_dt(declareRestartableData<Real>("unconstrained_dt", -1)),
    _at_sync_point(declareRestartableData<bool>("at_sync_point", false)),
    _first(declareRecoverableData<bool>("first", true)),
    _end_time(getParam<Real>("end_time")),
    _dtmin(getParam<Real>("dtmin")),
    _dtmax(getParam<Real>("dtmax")),
    _num_steps(getParam<unsigned int>("num_steps")),
    _n_startup_steps(getParam<int>("n_startup_steps")),
    _steps_taken(0),
    _trans_ss_check(getParam<bool>("trans_ss_check")),
    _ss_check_tol(getParam<Real>("ss_check_tol")),
    _ss_tmin(getParam<Real>("ss_tmin")),
    _old_time_solution_norm(declareRestartableData<Real>("old_time_solution_norm", 0.0)),
    _sync_times(_app.getOutputWarehouse().getSyncTimes()),
    _abort(getParam<bool>("abort_on_solve_fail")),
    _time_interval(declareRestartableData<bool>("time_interval", false)),
    _start_time(getParam<Real>("start_time")),
    _timestep_tolerance(getParam<Real>("timestep_tolerance")),
    _target_time(declareRestartableData<Real>("target_time", -1)),
    _use_multiapp_dt(getParam<bool>("use_multiapp_dt")),
    _verbose(getParam<bool>("verbose"))
{
  _problem.getNonlinearSystem().setDecomposition(_splitting);
  _t_step = 0;
  _dt = 0;
  _next_interval_output_time = 0.0;
  setSteps(_num_steps);

  // Either a start_time has been forced on us, or we want to tell the App about what our start time is (in case anyone else is interested.
  if (_app.hasStartTime())
    _start_time = _app.getStartTime();
  else
    _app.setStartTime(_start_time);

  _time = _time_old = _start_time;
  _problem.transient(true);

  if (parameters.isParamValid("predictor_scale"))
  {
    mooseWarning("Parameter 'predictor_scale' is deprecated, migrate your input file to use Predictor sub-block.");

    Real predscale = getParam<Real>("predictor_scale");
    if (predscale >= 0.0 && predscale <= 1.0)
    {
      InputParameters params = _app.getFactory().getValidParams("SimplePredictor");
      params.set<Real>("scale") = predscale;
      _problem.addPredictor("SimplePredictor", "predictor", params);
    }

    else
      mooseError("Input value for predictor_scale = "<< predscale << ", outside of permissible range (0 to 1)");
  }

  if (!_restart_file_base.empty())
    _problem.setRestartFile(_restart_file_base);

  setupTimeIntegrator();

  if (_app.halfTransient()) // Cut timesteps and end_time in half...
  {
    _end_time /= 2.0;
    _num_steps /= 2.0;

    if (_num_steps == 0) // Always do one step in the first half
      _num_steps = 1;
  }
}

Transient::~Transient()
{
}

void
Transient::init()
{
  if (!_time_stepper.get())
  {
    InputParameters pars = _app.getFactory().getValidParams("ConstantDT");
    pars.set<FEProblem *>("_fe_problem") = &_problem;
    pars.set<Transient *>("_executioner") = this;

    /**
     * We have a default "dt" set in the Transient parameters but it's possible for users to set other
     * parameters explicitly that could provide a better calculated "dt". Rather than provide difficult
     * to understand behavior using the default "dt" in this case, we'll calculate "dt" properly.
     */
    if (!_pars.isParamSetByAddParam("end_time") && !_pars.isParamSetByAddParam("num_steps") && _pars.isParamSetByAddParam("dt"))
      pars.set<Real>("dt") = (getParam<Real>("end_time") - getParam<Real>("start_time")) / static_cast<Real>(getParam<unsigned int>("num_steps"));
    else
      pars.set<Real>("dt") = getParam<Real>("dt");

    pars.set<bool>("reset_dt") = getParam<bool>("reset_dt");
    _time_stepper = MooseSharedNamespace::static_pointer_cast<TimeStepper>(_app.getFactory().create("ConstantDT", "TimeStepper", pars));
  }

  _problem.initialSetup();
  _time_stepper->init();

  if (_app.isRestarting())
    _time_old = _time;

  Moose::setup_perf_log.push("Output Initial Condition","Setup");
  _problem.outputStep(EXEC_INITIAL);
  Moose::setup_perf_log.pop("Output Initial Condition","Setup");

  // If this is the first step
  if (_t_step == 0)
    _t_step = 1;

  if (_t_step > 1) //Recover case
    _dt_old = _dt;

  else
  {
    computeDT();
    _dt = getDT();
  }


}

void
Transient::execute()
{
  std::cout<<"Junk"<<std::endl;

  takeStep();
}

void
Transient::computeDT()
{
  _time_stepper->computeStep(); // This is actually when DT gets computed
}

void
Transient::incrementStepOrReject()
{
  if (lastSolveConverged())
  {
    std::cout<<"Incrementing step!"<<std::endl;

#ifdef LIBMESH_ENABLE_AMR
    if (_problem.adaptivity().isOn())
      _problem.adaptMesh();
#endif

    _time_old = _time; // = _time_old + _dt;

    _problem.advanceState();
  }
  else
  {
    _problem.restoreMultiApps(EXEC_TIMESTEP_BEGIN);
    _problem.restoreMultiApps(EXEC_TIMESTEP_END);
    _time_stepper->rejectStep();
    _time = _time_old;
  }

  _first = false;

}

void
Transient::takeStep(Real input_dt)
{
  solveStep(input_dt);
}

void
Transient::solveStep(Real input_dt)
{
  _dt_old = _dt;

  if (input_dt == -1.0)
    _dt = computeConstrainedDT();
  else
    _dt = input_dt;

  Real current_dt = _dt;

  _problem.onTimestepBegin();

  // Increment time
  _time = _time_old + _dt;

  std::cout<<"Time: "<<_time<<std::endl;

  _time_stepper->preSolve();

  _problem.timestepSetup();


  _time_stepper->step();

  // We know whether or not the nonlinear solver thinks it converged, but we need to see if the executioner concurs
  if (lastSolveConverged())
  {
    _console << COLOR_GREEN << " Solve Converged!" << COLOR_DEFAULT << std::endl;

    if (getPicards() <= 1)
      _time_stepper->acceptStep();

    //output
    if (_time_interval && (_time + _timestep_tolerance >= _next_interval_output_time))
      _next_interval_output_time += _time_interval_output_interval;
  }
  else
  {
    _console << COLOR_RED << " Solve Did NOT Converge!" << COLOR_DEFAULT << std::endl;

    // Perform the output of the current, failed time step (this only occurs if desired)
    _fe_problem.outputStep(EXEC_FAILED);
  }
  _time_stepper->postSolve();
  _dt = current_dt; // _dt might be smaller than this at this point for multistep methods
}

void
Transient::beginStep()
{
  if (_first != true)
    incrementStepOrReject();

  _first = false;
}

void
Transient::endTransientStep(Real input_time)
{
  if (input_time == -1.0)
    _time = _time_old + _dt;
  else
    _time = input_time;
}

Real
Transient::computeConstrainedDT()
{
//  // If start up steps are needed
//  if (_t_step == 1 && _n_startup_steps > 1)
//    _dt = _input_dt/(double)(_n_startup_steps);
//  else if (_t_step == 1+_n_startup_steps && _n_startup_steps > 1)
//    _dt = _input_dt;

  Real dt_cur = _dt;
  std::ostringstream diag;

  //After startup steps, compute new dt
  if (_t_step > _n_startup_steps)
    dt_cur = getDT();

  else
  {
    diag << "Timestep < n_startup_steps, using old dt: "
         << std::setw(9)
         << std::setprecision(6)
         << std::setfill('0')
         << std::showpoint
         << std::left
         << _dt
         << " tstep: "
         << _t_step
         << " n_startup_steps: "
         << _n_startup_steps
         << std::endl;
  }
  _unconstrained_dt = dt_cur;

  if (_verbose)
    _console << diag.str();

  diag.str("");
  diag.clear();

  // Allow the time stepper to limit the time step
  _at_sync_point = _time_stepper->constrainStep(dt_cur);

  // Don't let time go beyond next time interval output if specified
  if ((_time_interval) &&
      (_time + dt_cur + _timestep_tolerance >= _next_interval_output_time))
  {
    dt_cur = _next_interval_output_time - _time;
    _at_sync_point = true;

    diag << "Limiting dt for time interval output at time: "
         << std::setw(9)
         << std::setprecision(6)
         << std::setfill('0')
         << std::showpoint
         << std::left
         << _next_interval_output_time
         << " dt: "
         << std::setw(9)
         << std::setprecision(6)
         << std::setfill('0')
         << std::showpoint
         << std::left
         << dt_cur
         << std::endl;
  }

  // Adjust to a target time if set
  if (_target_time > 0 && _time + dt_cur + _timestep_tolerance >= _target_time)
  {
    dt_cur = _target_time - _time;
    _at_sync_point = true;

    diag << "Limiting dt for target time: "
         << std::setw(9)
         << std::setprecision(6)
         << std::setfill('0')
         << std::showpoint
         << std::left
         << _next_interval_output_time
         << " dt: "
         << std::setw(9)
         << std::setprecision(6)
         << std::setfill('0')
         << std::showpoint
         << std::left
         << dt_cur
         << std::endl;
  }

  // Constrain by what the multi apps are doing
  Real multi_app_dt = _problem.computeMultiAppsDT(EXEC_TIMESTEP_BEGIN);
  if (_use_multiapp_dt || multi_app_dt < dt_cur)
  {
    dt_cur = multi_app_dt;
    _at_sync_point = false;
    diag << "Limiting dt for MultiApps: "
         << std::setw(9)
         << std::setprecision(6)
         << std::setfill('0')
         << std::showpoint
         << std::left
         << dt_cur
         << std::endl;
  }
  multi_app_dt = _problem.computeMultiAppsDT(EXEC_TIMESTEP_END);
  if (multi_app_dt < dt_cur)
  {
    dt_cur = multi_app_dt;
    _at_sync_point = false;
    diag << "Limiting dt for MultiApps: "
         << std::setw(9)
         << std::setprecision(6)
         << std::setfill('0')
         << std::showpoint
         << std::left
         << dt_cur
         << std::endl;
  }

  if (_verbose)
    _console << diag.str();

  return dt_cur;
}

Real
Transient::getDT()
{
  return _time_stepper->getCurrentDT();
}

bool
Transient::keepGoing()
{
  bool keep_going = !_problem.isSolveTerminationRequested();

  // Check for stop condition based upon steady-state check flag:
  if (lastSolveConverged() && _trans_ss_check == true && _time > _ss_tmin)
  {
    // Compute new time solution l2_norm
    Real new_time_solution_norm = _problem.getNonlinearSystem().currentSolution()->l2_norm();

    // Compute l2_norm relative error
    Real ss_relerr_norm = fabs(new_time_solution_norm - _old_time_solution_norm)/new_time_solution_norm;

    // Check current solution relative error norm against steady-state tolerance
    if (ss_relerr_norm < _ss_check_tol)
    {
      _console << "Steady-State Solution Achieved at time: " << _time << std::endl;
      //Output last solve if not output previously by forcing it
      keep_going = false;
    }
    else // Keep going
    {
      // Update solution norm for next time step
      _old_time_solution_norm = new_time_solution_norm;
      // Print steady-state relative error norm
      _console << "Steady-State Relative Differential Norm: " << ss_relerr_norm << std::endl;
    }
  }

  // Check for stop condition based upon number of simulation steps and/or solution end time:
  if (static_cast<unsigned int>(_t_step) > _num_steps)
    keep_going = false;

  if ((_time>_end_time) || (fabs(_time-_end_time)<=_timestep_tolerance))
    keep_going = false;

  if (!lastSolveConverged() && _abort)
  {
    _console << "Aborting as solve did not converge and input selected to abort" << std::endl;
    keep_going = false;
  }

  return keep_going;
}

void
Transient::estimateTimeError()
{
}

bool
Transient::lastSolveConverged()
{
  std::cout<<"LSC: "<<(multiAppsConverged() && _time_stepper->converged())<<std::endl;

  return multiAppsConverged() && _time_stepper->converged();
}

void
Transient::beginLoop()
{
  // Add time period start times to sync times
  const std::vector<TimePeriod *> time_periods = _problem.getTimePeriods();
  for (unsigned int i = 0; i < time_periods.size(); ++i)
    _time_stepper->addSyncTime(time_periods[i]->start());

  _time_stepper->preExecute();
}

void
Transient::endLoop()
{
  _time_stepper->postExecute();
}

Problem &
Transient::problem()
{
  return _problem;
}

void
Transient::setTargetTime(Real target_time)
{
  _target_time = target_time;
}

void
Transient::setupTimeIntegrator()
{
  if (_time_scheme.isValid() && _problem.hasTimeIntegrator())
    mooseError("You cannot specify time_scheme in the Executioner and independently add a TimeIntegrator to the system at the same time");

  if (!_problem.hasTimeIntegrator())
  {
    if (!_time_scheme.isValid())
      _time_scheme = "implicit-euler";

    // backwards compatibility
    std::string ti_str;

    switch (_time_scheme)
    {
    case 0: ti_str = "ImplicitEuler"; break;
    case 1: ti_str = "ExplicitEuler"; break;
    case 2: ti_str = "CrankNicolson"; break;
    case 3: ti_str = "BDF2"; break;
    case 4: ti_str = "RungeKutta2"; break;
    case 5: ti_str = "Dirk"; mooseError("Dirk requires parameters, please use the TimeIntegrator block instead of the \"scheme\" parameter."); break;
    default: mooseError("Unknown scheme"); break;
    }

    InputParameters params = _app.getFactory().getValidParams(ti_str);
    _problem.addTimeIntegrator(ti_str, ti_str, params);
  }
}

std::string
Transient::getTimeStepperName()
{
  if (_time_stepper.get())
    return demangle(typeid(*_time_stepper).name());
  else
    return std::string();
}
