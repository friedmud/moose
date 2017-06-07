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

#include "RayProblem.h"

// MOOSE Includes
#include "MooseMesh.h"
#include "DisplacedProblem.h"
#include "NonlinearSystem.h"
#include "AuxiliarySystem.h"

template <>
InputParameters
validParams<RayProblem>()
{
  InputParameters params = validParams<FEProblem>();
  params.addRequiredParam<unsigned int>("num_groups", "The number of energy groups to use");
  params.addParam<unsigned int>(
      "num_polar", 1, "The number of polar angles being used (only valid for a 2D problem)");

  params.addParam<unsigned int>(
      "max_rays", std::numeric_limits<unsigned int>::max(), "Maximum number of rays to run.");

  params.addParam<bool>("solve_ray", true, "Whether or not to solve the Ray problem");
  params.addParam<bool>("solve_fe", false, "Whether or not to actually solve the FE problem!");

  return params;
}

RayProblem::RayProblem(const InputParameters & params)
  : FEProblem(params),
    _num_groups(getParam<unsigned int>("num_groups")),
    _num_polar(getParam<unsigned int>("num_polar")),
    _solve_ray(getParam<bool>("solve_ray")),
    _solve_fe(getParam<bool>("solve_fe"))
{
  createRaySystem();

  // Determine the length of Rays needed to fully cross the domain
  _b_box = MeshTools::bounding_box(_mesh.getMesh());
  _domain_max_length =
      (_b_box.second - _b_box.first).norm() + 3 * TOLERANCE; // A little extra for good measure
}

RayProblem::~RayProblem() {}

void
RayProblem::initialSetup()
{
  FEProblem::initialSetup();

  _ray_system->initialSetup();
}

void
RayProblem::timestepSetup()
{
  FEProblem::timestepSetup();

  _ray_system->timestepSetup();
}

void
RayProblem::solve()
{
  _console << "0 Nonlinear |R|" << std::endl;

  /*
  for (unsigned int i = 0; i < _physics_iterations; i++)
  {
    if (_solve_ray)
    {
      // Solve the Ray problem
      _ray_system->solve();

      computeUserObjects(EXEC_CUSTOM, Moose::PRE_AUX);
      _aux->compute(EXEC_CUSTOM);
      computeUserObjects(EXEC_CUSTOM, Moose::POST_AUX);
    }

    // Solve the finite-element problem
    if (_solve_fe)
    {
      Real fe_norm = computeResidualL2Norm();

      _max_fe_norm = std::max(fe_norm, _max_fe_norm);

      _console << "\n FE Relative Residual Norm: " << fe_norm / _max_fe_norm << "\n" << std::endl;

      if (fe_norm / _max_fe_norm < _fe_active_tolerance)
      {
        _console << "Ray can enter Active" << std::endl;
        _fe_active = true;
      }
      else
      {
        _console << "Ray can't enter Active" << std::endl;
        _fe_active = false;
      }

      if (fe_norm / _max_fe_norm < _fe_tolerance)
      {
        _console << "FE Converged!" << std::endl;
        _fe_converged = true;
      }

      // Are we finished?
      if (// _fe_converged // && _ray_converged)
        break;

      FEProblem::solve();

      // Update the element average temperature
      computeAuxiliaryKernels(EXEC_TIMESTEP_END);
    }
  }
*/
}

bool
RayProblem::converged()
{
  return true;
}

void
RayProblem::subdomainSetup(SubdomainID subdomain, THREAD_ID tid)
{
  FEProblem::subdomainSetup(subdomain, tid);

  _ray_system->subdomainSetup(subdomain, tid);
}

void
RayProblem::RaySubdomainSetup(SubdomainID subdomain, THREAD_ID tid)
{
  _ray_system->subdomainSetup(subdomain, tid);
}

void
RayProblem::reinitElem(const Elem * elem, THREAD_ID tid)
{
  //  std::cout << "RayProblem::reinitElem()" << std::endl;

  _ray_system->reinitElem(elem, tid);

  FEProblem::reinitElem(elem, tid);
}

void
RayProblem::reinitElemPhys(const Elem * elem, std::vector<Point> phys_points_in_elem, THREAD_ID tid)
{
  _ray_system->subdomainSetup(elem->subdomain_id(), tid);

  FEProblem::reinitElemPhys(elem, phys_points_in_elem, tid);
}

void
RayProblem::RayReinitElem(const Elem * elem, THREAD_ID tid)
{
  _ray_system->reinitElem(elem, tid);
}

void
RayProblem::prepareMaterials(SubdomainID blk_id, THREAD_ID tid)
{
  _ray_system->subdomainSetup(blk_id, tid);

  FEProblem::prepareMaterials(blk_id, tid);
}

void
RayProblem::reinitMaterials(SubdomainID blk_id, THREAD_ID tid, bool swap_stateful)
{
  FEProblem::reinitMaterials(blk_id, tid, swap_stateful);
}

std::vector<VariableName>
RayProblem::getVariableNames()
{
  _console << "Getting variable names!" << std::endl;

  std::vector<VariableName> names = FEProblem::getVariableNames();

  const std::vector<VariableName> & ray_var_names = _ray_system->getVariableNames();
  names.insert(names.end(), ray_var_names.begin(), ray_var_names.end());

  return names;
}

bool
RayProblem::hasVariable(const std::string & var_name)
{
  return FEProblem::hasVariable(var_name) || _ray_system->hasVariable(var_name);
}

MooseVariable &
RayProblem::getVariable(THREAD_ID tid, const std::string & var_name)
{
  if (_ray_system->hasVariable(var_name))
    return _ray_system->getVariable(tid, var_name);

  return FEProblem::getVariable(tid, var_name);
}
