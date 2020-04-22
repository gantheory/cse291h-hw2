#ifndef _FLUID_H_
#define _FLUID_H_

#include <random>
#include "box.h"
#include "core.h"
#include "mesh.h"

class Fluid : public Mesh {
 private:
  // Constants for rendring small spheres.
  const float kRenderRadius = 0.04;
  const int stackCount = 5;
  const int sectorCount = 5;

  const glm::mat4 kModel = glm::mat4(1.0f);
  const glm::vec3 kColor = glm::vec3(1.0f);

  GLuint VAO;
  GLuint VBO_positions, VBO_normals, VBO_instance, EBO;

  float radius;
  std::vector<glm::vec3> positions;

  // Methods for rendering small spheres.
  std::vector<glm::vec3> surfaceVertices;
  std::vector<glm::vec3> surfaceNormals;
  std::vector<unsigned int> surfaceIndices;
  std::vector<glm::vec3> allSurfaceVertices;
  std::vector<glm::vec3> allSurfaceNormals;
  std::vector<unsigned int> allSurfaceIndices;

  glm::vec3 uniformSampleInsideSphere(float radius,
                                      std::default_random_engine& generator);

  void CreateSphereSurface();

  void AddPoints(glm::vec3 point);

 public:
  Fluid(Box& box, int numOfParticles);

  ~Fluid();

  void Draw(const glm::mat4& viewProjMtx, GLuint shader);

  void Update();

  glm::vec3 GetPosition(int i);

  int GetNumOfParticles();

  void SetPosition(int i, glm::vec3 p);

  float GetRadius();
};

#endif
