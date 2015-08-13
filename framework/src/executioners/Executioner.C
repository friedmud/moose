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

#include "Executioner.h"

// Moose includes
#include "MooseMesh.h"
#include "FEProblem.h"
#include "MooseApp.h"

// C++ includes
#include <vector>
#include <limits>

template<>
InputParameters validParams<Executioner>()
{
  InputParameters params = validParams<MooseObject>();
  params.addParam<FileNameNoExtension>("restart_file_base", "", "File base name used for restart");

  params.registerBase("Executioner");

  params.addParamNamesToGroup("restart_file_base", "Restart");

  params.addParam<std::vector<std::string> >("splitting", "Top-level splitting defining a hierarchical decomposition into subsystems to help the solver.");

  params.addParam<unsigned int>("picard_max_its", 1, "Number of times each timestep will be solved.  Mainly used when wanting to do Picard iterations with MultiApps that are set to execute_on picard_end, picard_begin, stage_end or stage_begin");
  params.addParam<Real>("picard_rel_tol", 1e-8, "The relative nonlinear residual drop to shoot for during Picard iterations.  This check is performed based on the Master app's nonlinear residual.");
  params.addParam<Real>("picard_abs_tol", 1e-50, "The absolute nonlinear residual to shoot for during Picard iterations.  This check is performed based on the Master app's nonlinear residual.");

#ifdef LIBMESH_HAVE_PETSC
  params += commonExecutionParameters();
#endif //LIBMESH_HAVE_PETSC

  params.addParamNamesToGroup("picard_max_its picard_rel_tol picard_abs_tol", "Picard");

  return params;
}

InputParameters commonExecutionParameters()
{
  InputParameters params = emptyInputParameters();
  MooseEnum solve_type("PJFNK JFNK NEWTON FD LINEAR");
  params.addParam<MooseEnum>   ("solve_type",      solve_type,
                                "PJFNK: Preconditioned Jacobian-Free Newton Krylov "
                                "JFNK: Jacobian-Free Newton Krylov "
                                "NEWTON: Full Newton Solve "
                                "FD: Use finite differences to compute Jacobian "
                                "LINEAR: Solving a linear problem");

  // Line Search Options
#ifdef LIBMESH_HAVE_PETSC
#if PETSC_VERSION_LESS_THAN(3,3,0)
  MooseEnum line_search("default cubic quadratic none basic basicnonorms", "default");
#else
  MooseEnum line_search("default shell none basic l2 bt cp", "default");
#endif
  std::string addtl_doc_str(" (Note: none = basic)");
#else
  MooseEnum line_search("default", "default");
  std::string addtl_doc_str("");
#endif
  params.addParam<MooseEnum>   ("line_search",     line_search, "Specifies the line search type" + addtl_doc_str);

#ifdef LIBMESH_HAVE_PETSC
  MultiMooseEnum common_petsc_options("", "", true);

  params.addParam<MultiMooseEnum>("petsc_options", common_petsc_options, "Singleton PETSc options");
  params.addParam<std::vector<std::string> >("petsc_options_iname", "Names of PETSc name/value pairs");
  params.addParam<std::vector<std::string> >("petsc_options_value", "Values of PETSc name/value pairs (must correspond with \"petsc_options_iname\"");
#endif //LIBMESH_HAVE_PETSC
  return params;
}

Executioner::Executioner(const InputParameters & parameters) :
    MooseObject(parameters),
    UserObjectInterface(parameters),
    PostprocessorInterface(parameters),
    Restartable(parameters, "Executioners"),
    _fe_problem(*parameters.getCheckedPointerParam<FEProblem *>("_fe_problem", "This might happen if you don't have a mesh")),
    _initial_residual_norm(std::numeric_limits<Real>::max()),
    _old_initial_residual_norm(std::numeric_limits<Real>::max()),
    _restart_file_base(getParam<FileNameNoExtension>("restart_file_base")),
    _splitting(getParam<std::vector<std::string> >("splitting")),
    _steps(1),
    _cycles(1),
    _picards(1),
    _stages(1),
    _current_step(_fe_problem._t_step),
    _current_cycle(_fe_problem._current_cycle),
    _current_picard(_fe_problem._current_picard),
    _current_stage(_fe_problem._current_stage),
    _total_executioner_loop_iterations(_fe_problem._total_executioner_loop_iterations),
    _picard_converged(declareRestartableData<bool>("picard_converged", false)),
    _picard_initial_norm(declareRestartableData<Real>("picard_initial_norm", 0.0)),
    _picard_timestep_begin_norm(declareRestartableData<Real>("picard_timestep_begin_norm", 0.0)),
    _picard_timestep_end_norm(declareRestartableData<Real>("picard_timestep_end_norm", 0.0)),
    _picard_rel_tol(getParam<Real>("picard_rel_tol")),
    _picard_abs_tol(getParam<Real>("picard_abs_tol")),
    _multiapps_converged(declareRestartableData<bool>("multiapps_converged", true)),
    _last_solve_converged(declareRestartableData<bool>("last_solve_converged", true))
{
  _current_step = 0;
}

Executioner::~Executioner()
{
}

void
Executioner::init()
{
  _fe_problem.initialSetup();

  Moose::setup_perf_log.push("Output Initial Condition","Setup");
  _fe_problem.outputStep(EXEC_INITIAL);
  Moose::setup_perf_log.pop("Output Initial Condition","Setup");

  // Call derived class init functions
  _init();
}

void
Executioner::justGo()
{
  // NOTE: if you remove this line, you will see a subset of tests failing. Those tests might have a wrong answer and might need to be regolded.
  // The reason is that we actually move the solution back in time before we actually start solving (which I think is wrong).  So this call here
  // is to maintain backward compatibility and so that MOOSE is giving the same answer.  However, we might remove this call and regold the test
  // in the future eventually.
  if (!_app.isRecovering())
    _fe_problem.advanceState();

  beforeSteps();
  executeSteps();
  afterSteps();
}

void
Executioner::executeSteps()
{
  std::cout<<"_steps+1 "<<_steps+1<<std::endl;

  for (_current_step = 1; _current_step < _steps + 1 && keepStepping(); _current_step++)
  {
    _total_executioner_loop_iterations++;

    std::cout<<"here"<<std::endl;

    beginStep();

    _fe_problem.computeUserObjects(EXEC_TIMESTEP_BEGIN, UserObjectWarehouse::PRE_AUX);

    _fe_problem.timestepSetup();

    _fe_problem.computeAuxiliaryKernels(EXEC_TIMESTEP_BEGIN);
    _fe_problem.computeUserObjects(EXEC_TIMESTEP_BEGIN, UserObjectWarehouse::POST_AUX);

    // Compute Pre-Aux User Objects (Timestep begin)
    _fe_problem.computeUserObjects(EXEC_TIMESTEP_BEGIN, UserObjectWarehouse::PRE_AUX);

    // Compute TimestepBegin AuxKernels
    _fe_problem.computeAuxiliaryKernels(EXEC_TIMESTEP_BEGIN);

    // Compute Post-Aux User Objects (Timestep begin)
    _fe_problem.computeUserObjects(EXEC_TIMESTEP_BEGIN, UserObjectWarehouse::POST_AUX);

    _fe_problem.execTransfers(EXEC_TIMESTEP_BEGIN);
    _multiapps_converged = _fe_problem.execMultiApps(EXEC_TIMESTEP_BEGIN, getPicards() == 1);

    if (!_multiapps_converged)
      return;

    // Perform output for timestep begin
    _fe_problem.outputStep(EXEC_TIMESTEP_BEGIN);

    executeCycles();

    _picard_converged=false;

    _last_solve_converged = lastSolveConverged();

    _fe_problem.onTimestepEnd();


    if (_last_solve_converged)
    {
      // Compute the Error Indicators and Markers
      _fe_problem.computeIndicatorsAndMarkers();

      // Output MultiApps if we were doing Picard iterations
      if (_picards > 1)
      {
        _fe_problem.advanceMultiApps(EXEC_TIMESTEP_BEGIN);
        _fe_problem.advanceMultiApps(EXEC_TIMESTEP_END);
      }

      _fe_problem.computeAuxiliaryKernels(EXEC_TIMESTEP_END);
      _fe_problem.computeUserObjects(EXEC_TIMESTEP_END, UserObjectWarehouse::POST_AUX);
      _fe_problem.execTransfers(EXEC_TIMESTEP_END);
      _multiapps_converged = _fe_problem.execMultiApps(EXEC_TIMESTEP_END, getPicards() == 1);

      if (!_multiapps_converged)
        return;

      // Perform the output of the current time step
      _fe_problem.outputStep(EXEC_TIMESTEP_END);

      _solution_change_norm = _fe_problem.solutionChangeNorm();

      _fe_problem.onTimestepEnd();
    }

    endStep();
  }
}

void
Executioner::executeCycles()
{
  std::cout<<"executeCycles cycles: "<<_cycles<<std::endl;

  for (_current_cycle = 0; _current_cycle < _cycles && keepCycling(); _current_cycle++)
  {
    if (_cycles > 1)
      _total_executioner_loop_iterations++;

    beginCycle();

    computeUserObjectsAndAuxiliaryKernels(EXEC_CYCLE_BEGIN);

    executePicards();

    _fe_problem.computeIndicatorsAndMarkers();

    if (_current_cycle + 1 != _cycles)
    {
      _fe_problem.outputStep(EXEC_TIMESTEP_END);
      _fe_problem.adaptMesh();
    }

    computeUserObjectsAndAuxiliaryKernels(EXEC_CYCLE_END);

    _fe_problem.outputStep(EXEC_CYCLE_END);

    endCycle();
  }
}

void
Executioner::executePicards()
{
  for (_current_picard = 0; _current_picard < _picards && keepPicarding(); _current_picard++)
  {
    if (_picards > 1)
      _total_executioner_loop_iterations++;

    beginPicard();

    computeUserObjectsAndAuxiliaryKernels(EXEC_PICARD_BEGIN);

    if (_picards > 1)
    {
      _fe_problem.backupMultiApps(EXEC_TIMESTEP_BEGIN);
      _fe_problem.backupMultiApps(EXEC_TIMESTEP_END);

      _console << "\nBeginning Picard Iteration " << _current_picard << "\n" << std::endl;

      Real current_norm = _fe_problem.computeResidualL2Norm();

      if (_current_picard == 0) // First Picard iteration - need to save off the initial nonlinear residual
      {
        _picard_initial_norm = current_norm;
        _console << "Initial Picard Norm: " << _picard_initial_norm << '\n';
      }
    }

    // For every iteration other than the first, we need to restore the state of the MultiApps
    if (_current_picard > 0)
    {
      _fe_problem.restoreMultiApps(EXEC_TIMESTEP_BEGIN);
      _fe_problem.restoreMultiApps(EXEC_TIMESTEP_END);
    }

    if (getPicards() > 1)
    {
      _picard_timestep_begin_norm = _fe_problem.computeResidualL2Norm();

      _console << "Picard Norm after TIMESTEP_BEGIN MultiApps: " << _picard_timestep_begin_norm << '\n';
    }

    executeStages();

    computeUserObjectsAndAuxiliaryKernels(EXEC_PICARD_END);


    if (_picards > 1)
    {
      _picard_timestep_end_norm = _fe_problem.computeResidualL2Norm();

      _console << "Picard Norm after TIMESTEP_END MultiApps: " << _picard_timestep_end_norm << '\n';

      Real max_norm = std::max(_picard_timestep_begin_norm, _picard_timestep_end_norm);

      Real max_relative_drop = max_norm / _picard_initial_norm;

      if (max_norm < _picard_abs_tol || max_relative_drop < _picard_rel_tol)
      {
        _console << "Picard converged!" << std::endl;

        _picard_converged = true;
        return;
      }
    }

    endPicard();
  }
}

void
Executioner::executeStages()
{
  for (_current_stage = 0; _current_stage < _stages && keepStaging(); _current_stage++)
  {
    if (_stages > 1)
      _total_executioner_loop_iterations++;

    beginStage();

    computeUserObjectsAndAuxiliaryKernels(EXEC_STAGE_BEGIN);

    execute();

    computeUserObjectsAndAuxiliaryKernels(EXEC_STAGE_END);

    endStage();
  }
}

void
Executioner::computeUserObjectsAndAuxiliaryKernels(ExecFlagType exec_flag)
{
  if (exec_flag == EXEC_CYCLE_END)
    std::cout<<"computeUserObjectsAndAuxiliaryKernels(CYCLE_END)"<<std::endl;

  _fe_problem.computeUserObjects(exec_flag, UserObjectWarehouse::PRE_AUX);

  _fe_problem.computeAuxiliaryKernels(exec_flag);
  _fe_problem.computeUserObjects(exec_flag, UserObjectWarehouse::POST_AUX);
}


std::string
Executioner::getTimeStepperName()
{
  return std::string();
}


void
Executioner::addAttributeReporter(const std::string & name, Real & attribute, const std::string execute_on)
{
  FEProblem * problem = parameters().getCheckedPointerParam<FEProblem *>("_fe_problem", "Failed to retrieve FEProblem when adding a attribute reporter in Executioner");
  InputParameters params = _app.getFactory().getValidParams("ExecutionerAttributeReporter");
  params.set<Real *>("value") = &attribute;
  if (!execute_on.empty())
    params.set<MultiMooseEnum>("execute_on") = execute_on;
  problem->addPostprocessor("ExecutionerAttributeReporter", name, params);
}

bool
Executioner::lastSolveConverged()
{
  return _last_solve_converged;
}

bool
Executioner::multiAppsConverged()
{
  return _multiapps_converged;
}
