#ifndef RAYTRACINGCOMMON_H
#define RAYTRACINGCOMMON_H

// libMesh Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/threads.h"

// System Includes
#include <chrono>

using namespace libMesh;

// tags
#define RAYS_STARTED 20000
#define RAYS_FINISHED 21000
#define RAY_BUFFER 22000
#define STOP_TRACING 23000

extern std::chrono::steady_clock::duration global_packing_time;
extern std::chrono::steady_clock::duration global_unpacking_time;

#endif // RAYTRACINGCOMMON_H
