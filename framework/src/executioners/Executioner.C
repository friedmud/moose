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

#ifdef LIBMESH_HAVE_PETSC
  params += commonExecutionParameters();
#endif //LIBMESH_HAVE_PETSC

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
    _current_step(_fe_problem.timeStep()),
    _current_cycle(declareRecoverableData<unsigned int>("current_cycle", 0)),
    _current_picard(declareRecoverableData<unsigned int>("current_picard", 0)),
    _current_stage(declareRecoverableData<unsigned int>("current_stage", 0)),
    _time(_fe_problem.time())
{
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
  for (_current_step = 1; _current_step < _steps + 1 && keepStepping(); _current_step++)
  {
    beginStep();

    switch (_time_scheme)
    {
      case REAL_TIME:
        _time += 0.1;
        break;
      case PSEUDO_TIME:
        _time = _current_step;
        break;
    }

    executeCycles();

    _fe_problem.outputStep(EXEC_TIMESTEP_END);

    endStep();
  }
}

void
Executioner::executeCycles()
{
  for (_current_cycle = 0; _current_cycle < _cycles && keepCycling(); _current_cycle++)
  {
    beginCycle();

    executePicards();

    _fe_problem.computeIndicatorsAndMarkers();

    if (_current_cycle + 1 != _cycles)
    {
      _fe_problem.outputStep(EXEC_TIMESTEP_END);
      _fe_problem.adaptMesh();
    }

    endCycle();
  }
}

void
Executioner::executePicards()
{
  for (_current_picard = 0; _current_picard < _picards && keepPicarding(); _current_picard++)
  {
    beginPicard();

    _fe_problem.computeUserObjects(EXEC_TIMESTEP_BEGIN, UserObjectWarehouse::PRE_AUX);

    _fe_problem.timestepSetup();
    _fe_problem.computeAuxiliaryKernels(EXEC_TIMESTEP_BEGIN);
    _fe_problem.computeUserObjects(EXEC_TIMESTEP_BEGIN, UserObjectWarehouse::POST_AUX);

    executeStages();

    _fe_problem.computeUserObjects(EXEC_TIMESTEP_END, UserObjectWarehouse::PRE_AUX);
    _fe_problem.onTimestepEnd();

    _fe_problem.computeAuxiliaryKernels(EXEC_TIMESTEP_END);
    _fe_problem.computeUserObjects(EXEC_TIMESTEP_END, UserObjectWarehouse::POST_AUX);

    endPicard();
  }
}

void
Executioner::executeStages()
{
  for (_current_stage = 0; _current_stage < _stages && keepStaging(); _current_stage++)
  {
    beginStage();
    execute();
    endStage();
  }
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
