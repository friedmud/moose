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

// Mocodile Includes
#include "PerProcessorRayTracingResultsVectorPostprocessor.h"
#include "PostprocessorInterface.h"
#include "RayProblem.h"

#include <algorithm>

registerMooseObject("MooseApp", PerProcessorRayTracingResultsVectorPostprocessor);

template <>
InputParameters
validParams<PerProcessorRayTracingResultsVectorPostprocessor>()
{
  InputParameters params = validParams<GeneralVectorPostprocessor>();

  MultiMooseEnum results(
      "rays_started rays_traced chunks_traced rays_received buffers_received "
      "rays_sent buffers_sent intersections generation_time propagation_time num_probes "
      "ray_pool_created receive_ray_pool_created receive_buffer_pool_created "
      "send_buffer_pool_created fast_lane_rays");

  params.addParam<MultiMooseEnum>("results", results, "The selection of results you want reported");

  return params;
}

PerProcessorRayTracingResultsVectorPostprocessor::PerProcessorRayTracingResultsVectorPostprocessor(
    const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _ray_problem(static_cast<RayProblem &>(_fe_problem)),
    _ray_tracing_study(_ray_problem.rayTracingStudy()),
    _results(getParam<MultiMooseEnum>("results")),
    _pid(processor_id()),
    _pid_values(declareVector("pid"))
{
  auto num_procs = n_processors();

  _pid_values.resize(n_processors(), 0);

  for (auto & result : _results)
  {
    std::string lower_result_name;

    std::transform(result.name().begin(),
                   result.name().end(),
                   std::back_inserter(lower_result_name),
                   ::tolower);

    auto id = result.id();

    _result_values[id] = &declareVector(lower_result_name);

    _result_values[id]->resize(num_procs, 0);
  }
}

void
PerProcessorRayTracingResultsVectorPostprocessor::initialize()
{
}

void
PerProcessorRayTracingResultsVectorPostprocessor::execute()
{
  _pid_values[_pid] = _pid;

  for (auto & result : _results)
  {
    switch (result)
    {
      case 0: // rays_started
        (*_result_values[0])[_pid] = _ray_tracing_study.localRaysStarted();
        break;
      case 1: // rays_traced
        (*_result_values[1])[_pid] = _ray_tracing_study.raysTraced();
        break;
      case 2: // chunks_traced
        (*_result_values[2])[_pid] = _ray_tracing_study.chunksTraced();
        break;
      case 3: // rays_received
        (*_result_values[3])[_pid] = _ray_tracing_study.raysReceived();
        break;
      case 4: // buffers_received
        (*_result_values[4])[_pid] = _ray_tracing_study.rayBuffersReceived();
        break;
      case 5: // rays_sent
        (*_result_values[5])[_pid] = _ray_tracing_study.raysSent();
        break;
      case 6: // buffers_sent
        (*_result_values[6])[_pid] = _ray_tracing_study.buffersSent();
        break;
      case 7: // intersections
        (*_result_values[7])[_pid] = _ray_tracing_study.localIntersections();
        break;
      case 8: // generation_time
        (*_result_values[8])[_pid] = _ray_tracing_study.generationTime();
        break;
      case 9: // propagation_time
        (*_result_values[9])[_pid] = _ray_tracing_study.propagationTime();
        break;
      case 10: // num_probes
        (*_result_values[10])[_pid] = _ray_tracing_study.numProbes();
        break;
      case 11: // ray_pool_created
        (*_result_values[11])[_pid] = _ray_problem._ray_pool.num_created();
        break;
      case 12: // receive_ray_pool_created
        (*_result_values[12])[_pid] = _ray_tracing_study.receiveRayPoolCreated();
        break;
      case 13: // receive_buffer_pool_created
        (*_result_values[13])[_pid] = _ray_tracing_study.receiveBufferPoolCreated();
        break;
      case 14: // send_buffer_pool_created
        (*_result_values[14])[_pid] = _ray_tracing_study.sendBufferPoolCreated();
        break;
      case 15: // fast_lane_rays
        (*_result_values[15])[_pid] = _ray_tracing_study.fastLaneRays();
        break;
      default:
        mooseError("Unknown result type '",
                   result.name(),
                   "/",
                   result.id(),
                   "' in PerProcessorRayTracingResultsVectorPostprocessor ",
                   name());
    }
  }
}

void
PerProcessorRayTracingResultsVectorPostprocessor::finalize()
{
  // TODO: Should be able to just "gather" to proc zero - but that's not working....
  gatherMax(_pid_values);

  for (auto & result : _results)
    gatherMax(*_result_values[result.id()]);
}
