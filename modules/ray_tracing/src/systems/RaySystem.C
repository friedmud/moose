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

#include "RaySystem.h"

// Local Includes
#include "RayProblem.h"
#include "RayKernel.h"
#include "RayBC.h"
#include "RayMaterial.h"

// MOOSE Includes
#include "MooseMesh.h"
#include "Factory.h"
#include "MooseTypes.h"

// libMesh Includes
#include "libmesh/quadrature_gauss.h"
#include "libmesh/node_range.h"
#include "libmesh/numeric_vector.h"

// System Includes
#include <numeric>

RaySystem::RaySystem(RayProblem & subproblem, const std::string & name, unsigned int num_groups)
  : SystemBase(subproblem, name, Moose::VAR_AUXILIARY),
    ConsoleStreamInterface(subproblem.getMooseApp()),
    _ray_problem(subproblem),
    _sys(_ray_problem.es().add_system<TransientExplicitSystem>(name)),
    _threaded_data(libMesh::n_threads()),
    _num_groups(num_groups),
    _current_group_solution(dynamic_cast<PetscVector<Number> &>(solution()))
{
  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
    _threaded_data[tid]._group_solution = dynamic_cast<PetscVector<Number> *>(
        &_sys.add_vector("thread_group_solution_" + std::to_string(tid), false, PARALLEL));

  FEType var_type(CONSTANT, MONOMIAL);

  for (unsigned int g = 0; g < _num_groups; g++)
    addVariable("group_" + std::to_string(g), var_type, 1.0);

  _console << "\nCreated RaySystem!!!\n" << std::endl;
}

RaySystem::~RaySystem() {}

void
RaySystem::init()
{
}

void
RaySystem::initialSetup()
{
  if (_ray_problem._solve_ray)
  {
    // Check for domain coverage for RayKernels
    {
      std::set<SubdomainID> subdomains;
      std::set<std::string> variables;
      _ray_kernels.subdomainsCovered(subdomains, variables);

      std::set<SubdomainID> missing;
      std::set_difference(_mesh.meshSubdomains().begin(),
                          _mesh.meshSubdomains().end(),
                          subdomains.begin(),
                          subdomains.end(),
                          std::inserter(missing, missing.begin()));

      if (!missing.empty())
      {
        std::ostringstream error;
        error << "Subdomains { ";
        std::copy(missing.begin(), missing.end(), std::ostream_iterator<SubdomainID>(error, " "));
        error << "} do not have RayKernels defined!";

        mooseError(error.str());
      }
    }

    // Check for domain coverage for RayMaterials
    {
      std::set<SubdomainID> subdomains;
      std::set<std::string> variables;
      _ray_materials.subdomainsCovered(subdomains, variables);

      std::set<SubdomainID> missing;
      std::set_difference(_mesh.meshSubdomains().begin(),
                          _mesh.meshSubdomains().end(),
                          subdomains.begin(),
                          subdomains.end(),
                          std::inserter(missing, missing.begin()));

      if (!missing.empty())
      {
        std::ostringstream error;
        error << "Subdomains { ";
        std::copy(missing.begin(), missing.end(), std::ostream_iterator<SubdomainID>(error, " "));
        error << "} do not have RayMaterials defined!";

        mooseError(error.str());
      }
    }
  }

  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
  {
    _ray_kernels.initialSetup(tid);
    _ray_materials.initialSetup(tid);
  }

  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
    updateActive(tid);
}

void
RaySystem::timestepSetup()
{
  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
  {
    _ray_kernels.timestepSetup(tid);
    _ray_materials.timestepSetup(tid);
  }
}

void
RaySystem::subdomainSetup()
{
  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
  {
    _ray_kernels.subdomainSetup(tid);
    _ray_materials.subdomainSetup(tid);
  }
}

void
RaySystem::jacobianSetup()
{
  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
  {
    _ray_kernels.jacobianSetup(tid);
    _ray_materials.jacobianSetup(tid);
  }
}

void
RaySystem::residualSetup()
{
  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
  {
    _ray_kernels.residualSetup(tid);
    _ray_materials.residualSetup(tid);
  }
}

void
RaySystem::updateActive(THREAD_ID tid)
{
  _ray_kernels.updateActive(tid);
  _ray_bcs.updateActive(tid);
  _ray_materials.updateActive(tid);
}

void
RaySystem::addRayKernel(const std::string & mk_name,
                        const std::string & name,
                        InputParameters parameters)
{
  _console << "Adding RayKernel!" << std::endl;

  parameters.set<FEProblem *>("_fe_problem") = &_ray_problem;
  parameters.set<FEProblemBase *>("_fe_problem_base") = &_ray_problem;
  parameters.set<SubProblem *>("_subproblem") = &_ray_problem;

  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); tid++)
  {
    MooseSharedPointer<RayKernel> mk = _factory.create<RayKernel>(mk_name, name, parameters, tid);
    _ray_kernels.addObject(mk, tid);
  }
}

void
RaySystem::addRayBC(const std::string & mbc_name,
                    const std::string & name,
                    InputParameters parameters)
{
  _console << "Adding RayBCl!" << std::endl;

  parameters.set<FEProblem *>("_fe_problem") = &_ray_problem;
  parameters.set<FEProblemBase *>("_fe_problem_base") = &_ray_problem;
  parameters.set<SubProblem *>("_subproblem") = &_ray_problem;

  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); tid++)
  {
    MooseSharedPointer<RayBC> mbc = _factory.create<RayBC>(mbc_name, name, parameters, tid);
    _ray_bcs.addObject(mbc, tid);
  }
}

void
RaySystem::addRayMaterial(const std::string & mk_name,
                          const std::string & name,
                          InputParameters parameters)
{
  _console << "Adding RayMaterial!" << std::endl;

  parameters.set<FEProblem *>("_fe_problem") = &_ray_problem;
  parameters.set<FEProblemBase *>("_fe_problem_base") = &_ray_problem;
  parameters.set<SubProblem *>("_subproblem") = &_ray_problem;

  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); tid++)
  {
    MooseSharedPointer<RayMaterial> mk =
        _factory.create<RayMaterial>(mk_name, name, parameters, tid);
    _ray_materials.addObject(mk, tid);
  }
}

void
RaySystem::subdomainSetup(SubdomainID current_subdomain, THREAD_ID tid)
{
  _threaded_data[tid]._current_ray_material =
      _ray_materials.getBlockObjects(current_subdomain, tid)[0];
}

void
RaySystem::augmentSparsity(SparsityPattern::Graph & /*sparsity*/,
                           std::vector<dof_id_type> & /*n_nz*/,
                           std::vector<dof_id_type> & /*n_oz*/)
{
}

void
RaySystem::solve()
{
  Moose::perf_log.push("RaySystem::solve()", "Execution");

  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
  {
    _ray_kernels.updateActive(tid);
    _ray_materials.updateActive(tid);
  }

  sweep();

  postSweep();

  Moose::perf_log.pop("RaySystem::solve()", "Execution");
}

void
RaySystem::reinitElem(const Elem * elem, THREAD_ID tid, bool only_sigma_t)
{
  auto & threaded_data = _threaded_data[tid];

  dof_id_type dof_number = elem->dof_number(number(), 0, 0);

  threaded_data._current_offset = _current_group_solution.map_global_to_local_index(dof_number);

  if (only_sigma_t)
    threaded_data._current_ray_material->reinitSigmaT(elem);
  else
    threaded_data._current_ray_material->reinit(elem);

  threaded_data._current_elem = elem;
}

void
RaySystem::sweep()
{
  VecGetArray(_current_group_solution.vec(), &_current_group_solution_values);
  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
    VecGetArray(_threaded_data[tid]._group_solution->vec(),
                &_threaded_data[tid]._group_solution_values);

  _ray_problem.rayTracingStudy().executeStudy();

  // Close up all of the vectors
  VecRestoreArray(_current_group_solution.vec(), &_current_group_solution_values);
  _current_group_solution.close();

  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
  {
    VecRestoreArray(_threaded_data[tid]._group_solution->vec(),
                    &_threaded_data[tid]._group_solution_values);

    _threaded_data[tid]._group_solution->close();
  }
}

void
RaySystem::postSweep()
{
  _current_group_solution.zero();

  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
    _current_group_solution += (*_threaded_data[tid]._group_solution);

  _current_group_solution.close();
}
