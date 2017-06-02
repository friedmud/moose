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
#include "RunStudy.h"
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

  params.addParam<unsigned int>("k_window_size", 20, "The size of the window for averaging k.");
  params.addParam<Real>(
      "k_window_tolerance",
      1e-3,
      "A change in the window value less than this will trigger the growth region.");
  params.addParam<unsigned int>(
      "k_window_bumps",
      3,
      "Number of times the k_window_tolerance must be met before entering the growth region.");
  params.addParam<Real>(
      "k_tolerance", 1e-5, "Terminate the solve whenever the active k changes by less than this.");

  params.addParam<Real>(
      "ray_growth",
      0.002,
      "Amaount to grow the number of rays by every iteration by during growth region.");
  params.addParam<Real>("ray_length_growth",
                        0.0002,
                        "Amaount to grow the ray length by every iteration during growth region.");

  params.addParam<Real>("growth_convergence_multiplier",
                        2.0,
                        "Sets the convergence criteria for moving from GROWTH to ACTIVE by "
                        "multiplying this number times the convergence criteria for kappa_source, "
                        "scalar_source and k.  Note: A higher number makes it 'easier' to go from "
                        "GROWTH to ACTIVE!");

  params.addParam<Real>(
      "scalar_flux_tolerance",
      5e-5,
      "Scalar flux relative norm must change by less than this before the solver quits.");

  params.addParam<Real>(
      "kappa_source_tolerance", 1e-4, "Convergence tolerance for kappa*sigma_f source");
  params.addParam<Real>("kappa_source_growth_tolerance",
                        1e-2,
                        "Convergence criteria for kappa*sigma_f to finish the growth region");

  params.addParam<bool>("solve_ray", true, "Whether or not to solve the Ray problem");
  params.addParam<bool>("solve_fe", false, "Whether or not to actually solve the FE problem!");

  params.addParam<bool>(
      "active_after_fe_converged",
      false,
      "Only allow entry into the active region after the other physics are converged");
  params.addParam<Real>("fe_active_tolerance",
                        3e-2,
                        "The active region can begin once the FE residual drops below this level");
  params.addParam<Real>("fe_tolerance", 1e-4, "Convergence tolerance for the FE problem");

  params.addParam<unsigned int>("power_iterations", 100, "Number of power iterations per step");

  params.addParam<unsigned int>(
      "physics_iterations", 1, "Number of (Picard) iterations per timestep");

  params.addParam<bool>(
      "terminate_after_ray_converged", true, "Stop the solve after the Ray problem has converged");

  params.addRequiredParam<UserObjectName>(
      "trrm_study", "The name of the RayTracingStudy UserObject to use for TRRM");

  params.addRequiredParam<UserObjectName>(
      "deterministic_study", "The name of the RayTracingStudy UserObject to use for Deterministic");

  params.addParam<bool>("trrm", true, "Whether or not we're doing TRRM");

  params.addParam<bool>(
      "volume_only", false, "Whether or not we're just doing a volume calculation");

  params.addParam<bool>("deterministic_active",
                        false,
                        "When true a TRRM solve will switch to a "
                        "deterministic solve once the active region "
                        "is entered");

  params.addParam<Real>(
      "dead_zone", 0, "How far a Ray should travel before modifying the scalar flux");

  params.addParamNamesToGroup(
      "num_groups num_polar max_rays k_window_size k_window_tolerance k_tolerance ray_growth "
      "ray_length_growth growth_convergence_multiplier scalar_flux_tolerance "
      "kappa_source_growth_tolerance solve_fe fe_tolerance fe_active_tolerance power_iterations "
      "physics_iterations active_after_fe_converged terminate_after_ray_converged trrm_study "
      "deterministic_active",
      "TRRM");

  params.addParamNamesToGroup("deterministic_study", "Deterministic");

  params.addParamNamesToGroup("num_groups num_polar power_iterations kappa_source_tolerance "
                              "physics_iterations trrm volume_only",
                              "Ray");

  params.addParam<unsigned int>("conditioning_iterations",
                                4,
                                "Number of power iterations to do (without updating the "
                                "source) when switching from TRRM to Deterministic");

  params.addPrivateParam<bool>("pke", false);

  params.addParamNamesToGroup("conditioning_iterations", "Hybrid");

  return params;
}

RayProblem::RayProblem(const InputParameters & params)
  : FEProblem(params),
    _num_groups(getParam<unsigned int>("num_groups")),
    _num_polar(getParam<unsigned int>("num_polar")),
    _trrm(getParam<bool>("trrm")),
    _deterministic_active(getParam<bool>("deterministic_active")),
    _volume_only(getParam<bool>("volume_only")),

    _max_rays(getParam<unsigned int>("max_rays")),

    _k_window_size(getParam<unsigned int>("k_window_size")),
    _k_window_tolerance(getParam<Real>("k_window_tolerance")),
    _k_window_bumps(getParam<unsigned int>("k_window_bumps")),
    _k_tolerance(getParam<Real>("k_tolerance")),

    _ray_growth(getParam<Real>("ray_growth")),
    _ray_length_growth(getParam<Real>("ray_length_growth")),

    _growth_convergence_multiplier(getParam<Real>("growth_convergence_multiplier")),

    _scalar_flux_tolerance(getParam<Real>("scalar_flux_tolerance")),

    _kappa_source_tolerance(getParam<Real>("kappa_source_tolerance")),
    _kappa_source_growth_tolerance(getParam<Real>("kappa_source_growth_tolerance")),

    _solve_ray(getParam<bool>("solve_ray")),

    _solve_fe(getParam<bool>("solve_fe")),
    _fe_tolerance(getParam<Real>("fe_tolerance")),
    _fe_active_tolerance(getParam<Real>("fe_active_tolerance")),

    _power_iterations(getParam<unsigned int>("power_iterations")),
    _conditioning_iterations(getParam<unsigned int>("conditioning_iterations")),

    _physics_iterations(getParam<unsigned int>("physics_iterations")),

    _active_after_converged(getParam<bool>("active_after_fe_converged")),
    _terminate_after_ray_converged(getParam<bool>("terminate_after_ray_converged")),

    _dead_zone(getParam<Real>("dead_zone")),

    //    _ray_system(std::make_shared<RaySystem>(*this, "Ray", _num_groups)),
    _ray_aux_system(*this, "RayAux", _num_groups),
    _pke(getParam<bool>("pke"))
{
  createRaySystem();

  // Determine the length of Rays needed to fully cross the domain
  _b_box = MeshTools::bounding_box(_mesh.getMesh());
  _domain_max_length =
      (_b_box.second - _b_box.first).norm() + 3 * TOLERANCE; // A little extra for good measure

  if (_mesh.dimension() == 3 && _num_polar != 1)
    mooseError("'Problem/num_polar' MUST be 1 for a 3D problem");
}

RayProblem::~RayProblem() {}

void
RayProblem::initialSetup()
{
  _current_ray_banks =
      std::make_shared<std::vector<std::vector<std::shared_ptr<Ray>>>>(libMesh::n_threads());
  _old_ray_banks =
      std::make_shared<std::vector<std::vector<std::shared_ptr<Ray>>>>(libMesh::n_threads());

  FEProblem::initialSetup();

  _ray_system->initialSetup();
}

void
RayProblem::timestepSetup()
{
  FEProblem::timestepSetup();

  _ray_converged = false;
  _fe_converged = false;
  _fe_active = false;

  _ray_system->timestepSetup();
}

void
RayProblem::solve()
{
  _console << "0 Nonlinear |R|" << std::endl;

  _max_fe_norm = 0;

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
      if (/*_fe_converged && */ _ray_converged)
        break;

      FEProblem::solve();

      // Update the element average temperature
      computeAuxiliaryKernels(EXEC_TIMESTEP_END);
    }
  }
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
  _ray_aux_system.subdomainSetup(subdomain, tid);
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
  _ray_aux_system.reinitElem(elem, tid);
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

  const std::vector<VariableName> & ray_aux_var_names = _ray_aux_system.getVariableNames();
  names.insert(names.end(), ray_aux_var_names.begin(), ray_aux_var_names.end());

  return names;
}

bool
RayProblem::hasVariable(const std::string & var_name)
{
  return FEProblem::hasVariable(var_name) || _ray_system->hasVariable(var_name) ||
         _ray_aux_system.hasVariable(var_name);
}

MooseVariable &
RayProblem::getVariable(THREAD_ID tid, const std::string & var_name)
{
  if (_ray_system->hasVariable(var_name))
    return _ray_system->getVariable(tid, var_name);

  if (_ray_aux_system.hasVariable(var_name))
    return _ray_aux_system.getVariable(tid, var_name);

  return FEProblem::getVariable(tid, var_name);
}
