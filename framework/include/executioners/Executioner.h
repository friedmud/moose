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

#ifndef EXECUTIONER_H
#define EXECUTIONER_H

#include "MooseObject.h"
#include "UserObjectInterface.h"
#include "PostprocessorInterface.h"
#include "Restartable.h"
#include "OutputWarehouse.h"

// System includes
#include <string>

class MooseMesh;
class Problem;
class Executioner;

template<>
InputParameters validParams<Executioner>();

/**
 * A helper function for creating execution related parameters, these are needed by
 * both Preconditioners and Executioners.
 */
InputParameters commonExecutionParameters();


/**
 * Executioners are objects that do the actual work of solving your problem.
 *
 * In general there are two "sets" of Executioners: Steady and Transient.
 *
 * The Difference is that a Steady Executioner usually only calls "solve()"
 * for the NonlinearSystem once... where Transient Executioners call solve()
 * multiple times... i.e. once per timestep.
 */
class Executioner :
  public MooseObject,
  public UserObjectInterface,
  public PostprocessorInterface,
  public Restartable
{
public:
  /**
   * Constructor
   *
   * @param name The name given to the Executioner in the input file.
   * @param parameters The parameters object holding data for the class to use.
   * @return Whether or not the solve was successful.
   */
  Executioner(const InputParameters & parameters);

  virtual ~Executioner();

  /**
   * Initialize the executioner
   *
   * DO NOT OVERRIDE THIS FUNCTION!
   *
   * Instead, override _init()
   */
  virtual void init();

  /**
   * "What are you preparing? You're always preparing!  Just Go!"
   *
   * The main exeuction loops.
   *
   * The execution loops are:
   *   for steps:
   *     for cycles:
   *       for picard:
   *         for stages:
   *           execute()
   *
   *   Each loop level has a before() and end() that can be overridden
   *
   * The main purpose of this function is to call executeSteps()
   */
  void justGo();

  /**
   * Get a reference to the FEProblem this Executioner is part of.
   */
  FEProblem & problem() { return _fe_problem; }

  /** The name of the TimeStepper
   * This is an empty string for non-Transient executioners
   * @return A string of giving the TimeStepper name
   */
  virtual std::string getTimeStepperName();

  /**
   * Can be used by subsclasses to call parentOutputPositionChanged()
   * on the underlying FEProblem.
   */
  virtual void parentOutputPositionChanged() {}

  /**
   * The current step.
   */
  int currentStep() { return _current_step; }

  /**
   * The current cycle.
   */
  int currentCycle() { return _current_cycle; }

  /**
   * The current picard.
   */
  int currentPicard() { return _current_picard; }

  /**
   * The current stage.
   */
  int currentStage() { return _current_stage; }

protected:
  /**
   * Called during init()
   *
   * Feel free to override this!
   */
  virtual void _init() {}


  /**
   * The primary method to override.  This is where the Executioner executes
   * the innermost loop
   */
  virtual void execute() = 0;

  /**
   * Set the number of steps to do
   */
  void setSteps(unsigned int steps) { _steps = steps; }

  /**
   * Set the number of cycles to do
   */
  void setCycles(unsigned int cycles) { _cycles = cycles; }

  /**
   * Set the number of picards to do
   */
  void setPicards(unsigned int picards) { _picards = picards; }

  /**
   * Set the number of stages to do
   */
  void setStages(unsigned int stages) { _stages = stages; }

  /**
   * Whether or not to continue the Steps loop
   */
  virtual bool keepStepping() { return true; }

  /**
   * Whether or not to continue the Cycles loop
   */
  virtual bool keepCycling() { return true; }

  /**
   * Whether or not to continue the Picards loop
   */
  virtual bool keepPicarding() { return true; }

  /**
   * Whether or not to continue the Stages loop
   */
  virtual bool keepStaging() { return true; }


  /**
   * Called at the beginning of each Step
   */
  virtual void beginStep() {}

  /**
   * Called at the beginning of each Cycle
   */
  virtual void beginCycle() {}

  /**
   * Called at the beginning of each Picard
   */
  virtual void beginPicard() {}

  /**
   * Called at the beginning of each Stage
   */
  virtual void beginStage() {}


  /**
   * Called at the ending of each Step
   */
  virtual void endStep() {}

  /**
   * Called at the ending of each Cycle
   */
  virtual void endCycle() {}

  /**
   * Called at the ending of each Picard
   */
  virtual void endPicard() {}

  /**
   * Called at the ending of each Stage
   */
  virtual void endStage() {}


  /**
   * Called before any looping starts
   */
  virtual void beforeSteps() {}

  /**
   * Called after all looping is finished
   */
  virtual void afterSteps() {}

  enum TimeScheme
  {
    REAL_TIME,
    PSEUDO_TIME
  };

  /**
   * Set the time computation scheme
   */
  void setTimeScheme(TimeScheme time_scheme) { _time_scheme = time_scheme; }

  /**
   * Adds a postprocessor to report a Real class attribute
   * @param name The name of the postprocessor to create
   * @param attribute The Real class attribute to report
   * @param execute_on When to execute the postprocessor that is created
   */
  virtual void addAttributeReporter(const std::string & name, Real & attribute, const std::string execute_on = "");

  FEProblem & _fe_problem;

  /// Initial Residual Variables
  Real _initial_residual_norm;
  Real _old_initial_residual_norm;

  // Restart
  std::string _restart_file_base;

  // Splitting
  std::vector<std::string> _splitting;

private:
  /**
   * The main solution steps
   * This loop will call executeCycles()
   */
  void executeSteps();

  /**
   * Cycles are mainly for mesh adaptivity loops
   * This loop will call executePicard()
   */
  void executeCycles();

  /**
   * Picard loops are for nonlinearly converging "tightly coupled" MultiApp simulations
   * This loop will call executeStages()
   */
  void executePicards();

  /**
   * Stages can be used for high order time integration or for Eigenvalue solves, etc.
   * This loop will call execute()
   */
  void executeStages();

  /**
   * Computes the UserObjects and AuxiliaryKernels for the given ExecFlagType
   */
  void computeUserObjectsAndAuxiliaryKernels(ExecFlagType exec_flag);

  /// The number of steps to do
  unsigned int _steps;

  /// The number of cycles to do
  unsigned int _cycles;

  /// The number of picard iterations to do
  unsigned int _picards;

  /// The number of stages to do
  unsigned int _stages;

  /// Current step
  int & _current_step;

  /// Current cycle
  unsigned int & _current_cycle;

  /// Current picard iteration
  unsigned int & _current_picard;

  /// Current stage
  unsigned int & _current_stage;

  /// The current time
  Real & _time;

  /// The current time computation scheme
  TimeScheme _time_scheme;
};

#endif //EXECUTIONER_H
