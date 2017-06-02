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
#include "RayBoundaryCondition.h"
#include "RayMaterial.h"
#include "RunStudy.h"

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
    _volume_only(subproblem.volumeOnly()),
    _fsr_sys(subproblem.es().add_system<ExplicitSystem>("fsr_system")),
    _fsr_sys_num(_fsr_sys.number()),
    _fsr_var(_fsr_sys.add_variable("fsr_var", FEType(CONSTANT, MONOMIAL))),
    _fsr_volumes(
        dynamic_cast<PetscVector<Number> &>(_fsr_sys.add_vector("fsr_volumes", false, PARALLEL))),
    _volume(computeVolume())
{
  for (unsigned int tid = 0; tid < libMesh::n_threads(); tid++)
    _threaded_data[tid]._fsr_volumes = dynamic_cast<PetscVector<Number> *>(
        &_fsr_sys.add_vector("thread_fsr_volumes_" + std::to_string(tid), false, PARALLEL));

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
RaySystem::addRayBoundaryCondition(const std::string & mbc_name,
                                   const std::string & name,
                                   InputParameters parameters)
{
  _console << "Adding RayBCl!" << std::endl;

  parameters.set<FEProblem *>("_fe_problem") = &_ray_problem;
  parameters.set<FEProblemBase *>("_fe_problem_base") = &_ray_problem;
  parameters.set<SubProblem *>("_subproblem") = &_ray_problem;

  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); tid++)
  {
    MooseSharedPointer<RayBoundaryCondition> mbc =
        _factory.create<RayBoundaryCondition>(mbc_name, name, parameters, tid);
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

  Moose::perf_log.pop("RaySystem::solve()", "Execution");
}

void
RaySystem::reinitElem(const Elem * elem, THREAD_ID tid, bool only_sigma_t)
{
  //  std::cout << "RaySystem::reinitElem()" << std::endl;

  auto & threaded_data = _threaded_data[tid];

  dof_id_type dof_number = elem->dof_number(number(), 0, 0);

  threaded_data._current_offset = _current_scalar_flux.map_global_to_local_index(dof_number);

  dof_id_type fsr_dof_number = elem->dof_number(_fsr_sys_num, 0, 0);

  threaded_data._current_fsr_offset = _fsr_volumes.map_global_to_local_index(fsr_dof_number);

  if (only_sigma_t)
    threaded_data._current_ray_material->reinitSigmaT(elem);
  else
    threaded_data._current_ray_material->reinit(elem);

  threaded_data._current_elem = elem;
}
