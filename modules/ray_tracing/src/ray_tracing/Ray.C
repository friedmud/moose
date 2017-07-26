#include "Ray.h"

bool
Ray::operator==(const Ray & other)
{
  return _data == other._data && _start == other._start && _end == other._end &&
         _starting_elem == other._starting_elem && _incoming_side == other._incoming_side &&
         _processor_crossings == other._processor_crossings &&
         _intersections == other._intersections && _distance == other._distance &&
         _integrated_distance == other._integrated_distance &&
         _azimuthal_spacing == other._azimuthal_spacing &&
         _azimuthal_weight == other._azimuthal_weight && _polar_spacing == other._polar_spacing &&
         _polar_sins == other._polar_sins && _polar_weights == other._polar_weights &&
         _ends_within_mesh == other._ends_within_mesh;
}
