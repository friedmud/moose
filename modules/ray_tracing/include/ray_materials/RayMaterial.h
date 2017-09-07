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
class RayProblemBase;

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
  RayProblemBase & _ray_problem;

  /// Number of energy groups
  unsigned int _num_groups;

  /// The thread ID this object is assigned to
  THREAD_ID _tid;
};

#endif
