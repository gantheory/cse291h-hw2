#ifndef _FLUIDDYNAMICS_H_
#define _FLUIDDYNAMICS_H_

#include "box.h"
#include "core.h"
#include "fluid.h"
#include "mesh.h"

class FluidDynamics {
 private:
  const float kDensity = 1.0;             // g / cm^3
  const float kGravity = 980.665f;        // cm / s^2
  const float kMaxDeltaT = 1e-4;          // s
  const float kKinematicViscosity = 1e0;  // cm^2 / s, 1e-0
  const float kStiffness = 1e-3;          // ???, 1e-3
  const float kSpringWall = 1e5;

  Fluid* fluid;
  Box* box;

  std::vector<float> mass;            // g
  std::vector<glm::vec3> velocities;  // cm / s
  std::vector<float> densities;       // g / cm^3

  int iteration;
  float deltaT;  // s
  float h;       // cm^3

  float W(int i, int j);

  glm::vec3 ComputeNablaW(int i, int j);

 public:
  FluidDynamics(Fluid* fluid, Box* box);

  void Update();
};

#endif
