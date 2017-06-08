#ifndef ELEMENTCENTROIDRAYS_H
#define ELEMENTCENTROIDRAYS_H

#include "RayTracingStudy.h"

class ElementCentroidRays;

template <>
InputParameters validParams<ElementCentroidRays>();

class ElementCentroidRays : public RayTracingStudy
{
public:
  ElementCentroidRays(const InputParameters & parameters);

protected:
  virtual void generateRays();

  const std::vector<dof_id_type> & _start_elems;
  const std::vector<dof_id_type> & _end_elems;
};

#endif
