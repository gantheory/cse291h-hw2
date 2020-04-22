#ifndef _FLUIDDYNAMICS_H_
#define _FLUIDDYNAMICS_H_

#include <unordered_map>
#include "box.h"
#include "core.h"
#include "fluid.h"
#include "mesh.h"

class FluidDynamics {
 private:
  const float kDensity = 1.0;              // g / cm^3
  const float kGravity = 980.665f;         // cm / s^2
  const float kMaxDeltaT = 1e-4;           // s
  const float kKinematicViscosity = 0.50;  // cm^2 / s, 1
  const float kStiffness = 5e-1;           // ???, 5
  const float kSpringWall = 1e5;
  const float h = 0.065;  // cm^3

  const float kLambda = 0.006;
  const int p1 = 73856093;
  const int p2 = 19349663;
  const int p3 = 83492791;
  const int kUpdateIter = 2;

  Fluid* fluid;
  Box* box;

  std::vector<float> mass;            // g
  std::vector<glm::vec3> velocities;  // cm / s
  std::vector<float> densities;       // g / cm^3

  int iteration;
  float deltaT;  // s

  std::vector<std::vector<int>> neighbors;
  std::unordered_map<int, std::vector<int>> spatialHash;

  float W(int i, int j);

  glm::vec3 ComputeNablaW(int i, int j);

  std::vector<int> GetNeighborsHashValue(glm::vec3& point);

  int hashValue(glm::vec3& point);

  void UpdateSpatialHashTable();

 public:
  FluidDynamics(Fluid* fluid, Box* box);

  void Update();
};

#endif
