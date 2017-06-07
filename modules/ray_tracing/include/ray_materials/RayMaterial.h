#ifndef RAYMATERIAL_H
#define RAYMATERIAL_H

// MOOSE Includes
#include "MooseObject.h"
#include "SetupInterface.h"
#include "Restartable.h"
#include "BlockRestrictable.h"
#include "Coupleable.h"
#include "UserObjectInterface.h"

// System Includes
#include <vector>

class RayMaterial;
class RayProblem;

template <>
InputParameters validParams<RayMaterial>();

/**
 * Holds all of the XS data for one material
 */
class RayMaterial : public MooseObject,
                    public SetupInterface,
                    public Restartable,
                    public BlockRestrictable,
                    public Coupleable,
                    public UserObjectInterface,
                    public TransientInterface
{
public:
  RayMaterial(const InputParameters & params);

  /// Called on each segment so the material can recompute the sigma_t
  virtual void reinitSigmaT(const Elem * /*elem*/){};

  /// Called on each segment so the material can recompute all XS
  virtual void reinit(const Elem * /*elem*/){};

  /// Whether or not this material produces fissions
  virtual bool fissionable(const Elem * /*elem*/) { return _fissionable; }

  /// Number of energy groups
  unsigned int _n_groups;

  /// Sigma total for each energy group (starting with highest energy)
  std::vector<double> _sigma_t;

  /// Sigma-Fission for each energy group (starting with highest energy)
  std::vector<double> _sigma_f;

  /// Nu-Sigma-Fission for each energy group (starting with highest energy)
  std::vector<double> _nu_sigma_f;

  /// Kappa-Sigma-Fission for each energy group (starting with highest energy)
  std::vector<double> _kappa_sigma_f;

  /// Fission spectrum for each energy group (starting with highest energy)
  std::vector<double> _chi;

  /// Scattering.  Out is the first index.  In is the second
  /// So _sigma_s[2][1] is the scattering cross section for group1 to group2
  /// It's done this way to make it simple to compute the group's contribution through a matvec
  std::vector<std::vector<double>> _sigma_s;

  /**
   * This needs to be here to satisfy the MooseObjectWarehouse interface
   */
  class DummyVariable
  {
  public:
    std::string name() { return ""; }
  };

  /**
   * This needs to be here to satisfy the MooseObjectWarehouse interface
   */
  DummyVariable variable() { return DummyVariable(); }

protected:
  /// Whether or not this material produces fissions
  bool _fissionable;

  RayProblem & _ray_problem;

  /// The thread ID this object is assigned to
  THREAD_ID _tid;
};

#endif
