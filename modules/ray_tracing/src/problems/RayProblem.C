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
  params.addRequiredParam<unsigned int>("num_groups", "The number of value spots on a Ray");
  params.addParam<unsigned int>(
      "num_polar", 1, "The number of polar angles being used (only valid for a 2D problem)");

  params.addParam<unsigned int>(
      "max_rays", std::numeric_limits<unsigned int>::max(), "Maximum number of rays to run.");

  params.addParam<bool>("solve_ray", true, "Whether or not to solve the Ray problem");
  params.addParam<bool>("solve_fe", false, "Whether or not to actually solve the FE problem!");

  params.addParam<UserObjectName>("study", "The RayTracing study");

  params.addParam<bool>(
      "ray_kernel_coverage_check", true, "Whether or not to check for coverage of RayKernels");

  params.addParam<bool>(
      "ray_material_coverage_check", true, "Whether or not to check for coverage of RayMaterials");

  return params;
}

RayProblemBase::RayProblemBase(const InputParameters & params)
  : FEProblem(params),
    _num_groups(getParam<unsigned int>("num_groups")),
    _num_polar(getParam<unsigned int>("num_polar")),
    _solve_ray(getParam<bool>("solve_ray")),
    _solve_fe(getParam<bool>("solve_fe")),
    _ray_kernel_coverage_check(getParam<bool>("ray_kernel_coverage_check")),
    _ray_material_coverage_check(getParam<bool>("ray_material_coverage_check"))
{
  // Determine the length of Rays needed to fully cross the domain
  _b_box = MeshTools::bounding_box(_mesh.getMesh());
  _domain_max_length =
      (_b_box.second - _b_box.first).norm() + 3 * TOLERANCE; // A little extra for good measure
}

RayProblemBase::~RayProblemBase() {}

void
RayProblemBase::initialSetup()
{
  FEProblem::initialSetup();

  _ray_system->initialSetup();
}

void
RayProblemBase::timestepSetup()
{
  FEProblem::timestepSetup();

  _ray_system->timestepSetup();
}

void
RayProblemBase::solve()
{
  _console << "0 Nonlinear |R|" << std::endl;

  // Solve the finite-element problem
  if (_solve_fe)
  {
    FEProblem::solve();
  }

  if (_solve_ray)
  {
    // Solve the Ray problem
    _ray_system->solve();
  }
}

bool
RayProblemBase::converged()
{
  return true;
}

void
RayProblemBase::subdomainSetup(SubdomainID subdomain, THREAD_ID tid)
{
  FEProblem::subdomainSetup(subdomain, tid);

  _ray_system->subdomainSetup(subdomain, tid);
}

void
RayProblemBase::RaySubdomainSetup(SubdomainID subdomain, THREAD_ID tid)
{
  _ray_system->subdomainSetup(subdomain, tid);
}

void
RayProblemBase::reinitElem(const Elem * elem, THREAD_ID tid)
{
  _ray_system->reinitElem(elem, tid);

  FEProblem::reinitElem(elem, tid);
}

void
RayProblemBase::reinitElemPhys(const Elem * elem,
                               std::vector<Point> phys_points_in_elem,
                               THREAD_ID tid)
{
  _ray_system->subdomainSetup(elem->subdomain_id(), tid);

  FEProblem::reinitElemPhys(elem, phys_points_in_elem, tid);
}

void
RayProblemBase::reinitNeighbor(const Elem * elem, unsigned int side, THREAD_ID tid)
{
  const Elem * neighbor = elem->neighbor_ptr(side);
  unsigned int neighbor_side = neighbor->which_neighbor_am_i(elem);

  _ray_system->reinitNeighborFace(neighbor, neighbor_side, 0, tid);

  FEProblem::reinitNeighbor(elem, side, tid);
}

void
RayProblemBase::RayReinitElem(const Elem * elem, THREAD_ID tid)
{
  _ray_system->reinitElem(elem, tid);
}

void
RayProblemBase::prepareMaterials(SubdomainID blk_id, THREAD_ID tid)
{
  _ray_system->subdomainSetup(blk_id, tid);

  FEProblem::prepareMaterials(blk_id, tid);
}

void
RayProblemBase::reinitMaterials(SubdomainID blk_id, THREAD_ID tid, bool swap_stateful)
{
  FEProblem::reinitMaterials(blk_id, tid, swap_stateful);
}

std::vector<VariableName>
RayProblemBase::getVariableNames()
{
  _console << "Getting variable names!" << std::endl;

  std::vector<VariableName> names = FEProblem::getVariableNames();

  const std::vector<VariableName> & ray_var_names = _ray_system->getVariableNames();
  names.insert(names.end(), ray_var_names.begin(), ray_var_names.end());

  return names;
}

bool
RayProblemBase::hasVariable(const std::string & var_name)
{
  return FEProblem::hasVariable(var_name) || _ray_system->hasVariable(var_name);
}

MooseVariable &
RayProblemBase::getVariable(THREAD_ID tid, const std::string & var_name)
{
  if (_ray_system->hasVariable(var_name))
    return _ray_system->getVariable(tid, var_name);

  return FEProblem::getVariable(tid, var_name);
}

RayProblem::RayProblem(const InputParameters & params) : RayProblemBase(params)
{
  _ray_system = std::make_shared<RaySystem>(*this, "Ray", _num_groups);
}
