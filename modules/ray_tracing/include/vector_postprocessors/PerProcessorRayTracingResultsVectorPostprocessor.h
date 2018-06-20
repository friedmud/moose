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

#ifndef PERPROCESSORRAYTRACINGRESULTSVECTORPOSTPROCESSOR_H
#define PERPROCESSORRAYTRACINGRESULTSVECTORPOSTPROCESSOR_H

#include "GeneralVectorPostprocessor.h"

// Forward Declarations
class RayProblem;
class RayTracingStudy;
class PerProcessorRayTracingResultsVectorPostprocessor;

template <>
InputParameters validParams<PerProcessorRayTracingResultsVectorPostprocessor>();

/**
 * Outputs per-processor metrics from a RayTracingStudy
 */
class PerProcessorRayTracingResultsVectorPostprocessor : public GeneralVectorPostprocessor
{
public:
  PerProcessorRayTracingResultsVectorPostprocessor(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;

protected:
  /// The MOCSystem to pull the data from
  RayProblem & _ray_problem;

  /// The study to pull data from
  RayTracingStudy & _ray_tracing_study;

  /// The results that were chosen
  const MultiMooseEnum & _results;

  /// For convenience and speed
  processor_id_type _pid;

  /// Vector to hold PIDs
  VectorPostprocessorValue & _pid_values;

  /// The VPP values
  std::map<int, VectorPostprocessorValue *> _result_values;
};

#endif
