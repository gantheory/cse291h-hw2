#include "fluid_dynamics.h"
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <random>
#include "core.h"

#define PI 3.1415926f
#define DEBUG_ITER 500
#define LAMBDA 0.006

FluidDynamics::FluidDynamics(Fluid* fluid, Box* box) : fluid(fluid), box(box) {
  int numOfParticles = fluid->GetNumOfParticles();

  // h = pow(4.f * fluid->GetRadius() / (3.f * (float)numOfParticles), 1.f
  // / 3.f);
  h = 0.15;

  // mass.resize(numOfParticles, PI * pow(h, 3) * kDensity);
  mass.resize(numOfParticles, pow(h, 3) * kDensity);

  velocities.resize(numOfParticles, glm::vec3(0.f));

  densities.resize(numOfParticles, kDensity);

  iteration = 0;
}

void checkNan(glm::vec3 p) {
  if (isnan(p.x) || isnan(p.y) || isnan(p.z)) {
    std::cerr << "nan vector found: " << to_string(p) << std::endl;
    assert(false);
  }
};

void FluidDynamics::Update() {
  ++iteration;
  if (iteration % DEBUG_ITER == 1) {
    std::cerr << "============================================================="
                 "==================="
              << std::endl;
  }

  int numOfParticles = fluid->GetNumOfParticles();
  glm::vec3 boxMin = box->GetBoxMin();
  glm::vec3 boxMax = box->GetBoxMax();
  std::vector<glm::vec3> viscosityForce(numOfParticles, glm::vec3(0.f));
  std::vector<glm::vec3> gravity(numOfParticles, glm::vec3(0.f));
  std::vector<glm::vec3> pressureForce(numOfParticles, glm::vec3(0.f));

  // Determine time step.
  float maxVelocity = -1e9;
  for (int i = 0; i < numOfParticles; ++i) {
    maxVelocity = fmax(maxVelocity, velocities[i].x);
    maxVelocity = fmax(maxVelocity, velocities[i].y);
    maxVelocity = fmax(maxVelocity, velocities[i].z);
  }
  deltaT = fmin(kMaxDeltaT, LAMBDA * h / maxVelocity);

  if (iteration % DEBUG_ITER == 1) {
    std::cerr << "Current time step (deltaT): " << deltaT << std::endl;
  }

  // Find neighbors.
  int total = 0;
  std::vector<std::vector<int>> neighbors(numOfParticles);
  for (int i = 0; i < numOfParticles; ++i) {
    glm::vec3 currentPoint = fluid->GetPosition(i);
    for (int j = 0; j < numOfParticles; ++j) {
      glm::vec3 anotherPoint = fluid->GetPosition(j);
      if (glm::distance(currentPoint, anotherPoint) < 2.f * h) {
        ++total;
        neighbors[i].push_back(j);
      }
    }
  }
  if (iteration % DEBUG_ITER == 1) {
    std::cerr << "Average number of neighbors: "
              << (float)total / (float)numOfParticles << std::endl;
  }

  // Compute nabla W.
  std::vector<std::vector<glm::vec3>> nablaW(numOfParticles);
  for (int i = 0; i < numOfParticles; ++i) {
    nablaW[i].resize(neighbors[i].size());
    for (int j = 0; j < (int)nablaW[i].size(); ++j)
      nablaW[i][j] = ComputeNablaW(i, neighbors[i][j]);
  }

  // viscosity
  for (int i = 0; i < numOfParticles; ++i) {
    glm::vec3 xi = fluid->GetPosition(i);
    glm::vec3 sum(0.f);
    for (int ii = 0; ii < (int)neighbors[i].size(); ++ii) {
      int j = neighbors[i][ii];
      glm::vec3 xj = fluid->GetPosition(j);
      glm::vec3 xij = xi - xj;
      sum += mass[j] / densities[j] * (velocities[i] - velocities[j]) *
             dot(xij, nablaW[i][ii]) /
             (float)(dot(xij, xij) + 0.01 * pow(h, 2));
    }
    viscosityForce[i] += mass[i] * kKinematicViscosity * 2 * sum;

    auto& p = viscosityForce[i];
    if (isnan(p.x) || isnan(p.y) || isnan(p.z)) {
      std::cerr << "nan viscosity force found: " << to_string(p) << std::endl;
      assert(false);
    }
  }

  // Gravity
  for (int i = 0; i < numOfParticles; ++i) gravity[i].y -= mass[i] * kGravity;

  // Splitting
  std::vector<float> pressure(numOfParticles);
  for (int i = 0; i < numOfParticles; ++i) {
    velocities[i] += deltaT * (viscosityForce[i] + gravity[i]) / mass[i];

    auto& p = velocities[i];
    if (isnan(p.x) || isnan(p.y) || isnan(p.z)) {
      std::cerr << "nan velocity found: " << to_string(p) << std::endl;
      assert(false);
    }
  }
  float averageDensity = 0.0;
  for (int i = 0; i < numOfParticles; ++i) {
    densities[i] = 0;
    for (int ii = 0; ii < (int)neighbors[i].size(); ++ii) {
      int j = neighbors[i][ii];
      densities[i] +=
          mass[j] * W(i, j) +
          deltaT * mass[j] * dot(velocities[i] - velocities[j], nablaW[i][ii]);
    }

    averageDensity += densities[i];
  }
  if (iteration % DEBUG_ITER == 1) {
    std::cerr << "Average density: " << averageDensity / (float)numOfParticles
              << std::endl;
  }

  // Pressure
  for (int i = 0; i < numOfParticles; ++i) {
    pressure[i] = kStiffness * (pow(densities[i] / kDensity, 7) - 1.f);
    if (isnan(pressure[i])) {
      std::cerr << "nan pressure found: " << pressure[i] << std::endl;
    }
  }
  for (int i = 0; i < numOfParticles; ++i) {
    glm::vec3 nablaP(0.f);
    for (int ii = 0; ii < (int)neighbors[i].size(); ++ii) {
      int j = neighbors[i][ii];
      nablaP += mass[j] *
                (float)(pressure[i] / pow(densities[i], 2) +
                        pressure[j] / pow(densities[j], 2)) *
                nablaW[i][ii];
    }
    pressureForce[i] += -mass[i] * nablaP;

    auto& p = pressureForce[i];
    if (isnan(p.x) || isnan(p.y) || isnan(p.z)) {
      std::cerr << "nan pressure force found: " << to_string(p) << std::endl;
      std::cerr << "pressure: " << pressure[i] << std::endl;
      std::cerr << "densities: " << densities[i] << std::endl;
      assert(false);
    }
  }

  // Collision to the ground and walls
  for (int i = 0; i < numOfParticles; ++i) {
    glm::vec3 currentPoint = fluid->GetPosition(i);
    if (currentPoint.x < boxMin.x) {
      float force = kSpringWall * (boxMin.x - currentPoint.x);
      velocities[i].x += deltaT * force / mass[i];
    } else if (currentPoint.x > boxMax.x) {
      float force = kSpringWall * (currentPoint.x - boxMax.x);
      velocities[i].x -= deltaT * force / mass[i];
    }
    if (currentPoint.y < boxMin.y) {
      float force = kSpringWall * (boxMin.y - currentPoint.y);
      velocities[i].y += deltaT * force / mass[i];
    } else if (currentPoint.y > boxMax.y) {
      float force = kSpringWall * (currentPoint.y - boxMax.y);
      velocities[i].y -= deltaT * force / mass[i];
    }
    if (currentPoint.z < boxMin.z) {
      float force = kSpringWall * (boxMin.z - currentPoint.z);
      velocities[i].z += deltaT * force / mass[i];
    } else if (currentPoint.z > boxMax.z) {
      float force = kSpringWall * (currentPoint.z - boxMax.z);
      velocities[i].z -= deltaT * force / mass[i];
    }
  }

  // Integration.
  for (int i = 0; i < numOfParticles; ++i) {
    glm::vec3 acceleration = pressureForce[i] / mass[i];
    velocities[i] += acceleration * deltaT;
    glm::vec3 currentPoint = fluid->GetPosition(i);
    currentPoint += velocities[i] * deltaT;
    if (isnan(currentPoint.x) || isnan(currentPoint.y) ||
        isnan(currentPoint.z)) {
      std::cerr << "nan position found: " << to_string(currentPoint)
                << std::endl;
      assert(false);
    }
    fluid->SetPosition(i, currentPoint);
  }

  if (iteration % DEBUG_ITER == 1) {
    std::cerr << "============================================================="
                 "==================="
              << std::endl;
  }
}

float FluidDynamics::W(int i, int j) {
  glm::vec3 xi = fluid->GetPosition(i);
  glm::vec3 xj = fluid->GetPosition(j);

  float q = glm::distance(xi, xj) / h;

  float fq = 3.f / (2.f * PI);
  if (0 <= q && q < 1) {
    fq *= (2.f / 3.f - pow(q, 2) + 0.5f * pow(q, 3));
  } else if (1 <= q && q < 2) {
    fq *= (1.f / 6.f * pow(2 - q, 3));
  } else {
    fq = 0.f;
  }
  return fq / pow(h, 3);
}

glm::vec3 FluidDynamics::ComputeNablaW(int i, int j) {
  glm::vec3 xi = fluid->GetPosition(i);
  glm::vec3 xj = fluid->GetPosition(j);

  float q = glm::distance(xi, xj) / h;

  float k = 3.f / (2.f * PI * pow(h, 3));
  if (0 <= q && q < 1) {
    return k *
           (float)(-2.f / pow(h, 2) +
                   3.f / 2.f * glm::distance(xi, xj) / pow(h, 3)) *
           (xi - xj);
  } else if (1 <= q && q < 2) {
    float WPartialQ = -0.5f * pow(2.f - q, 2);
    glm::vec3 qPartialXi = (xi - xj) / (h * glm::distance(xi, xj));
    return k * WPartialQ * qPartialXi;
  } else if (2 <= q) {
    return glm::vec3(0.f);
  }
  assert(false);
}
