#include "fluid.h"
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <random>
#include "core.h"

#define PI 3.1415926f

Fluid::Fluid(Box& box, int numOfParticles) {
  // Sample points
  glm::vec3 boxMin = box.GetBoxMin();
  glm::vec3 boxMax = box.GetBoxMax();
  radius =
      fmin((boxMax.x - boxMin.x) / 2 * 0.95, (boxMax.z - boxMin.z) / 2 * 0.95);
  glm::vec3 initialCenter((boxMin.x + boxMax.x) / 2, boxMax.y - radius,
                          (boxMin.z + boxMax.z) / 2);
  std::cerr << "Initial radius: " << radius << std::endl;
  std::cerr << "Initial center: " << to_string(initialCenter) << std::endl;

  CreateSphereSurface();
  std::default_random_engine generator;
  for (int i = 0; i < numOfParticles; ++i) {
    glm::vec3 center = uniformSampleInsideSphere(radius, generator);
    center += initialCenter;
    positions.push_back(center);

    allSurfaceVertices.insert(allSurfaceVertices.end(), surfaceVertices.begin(),
                              surfaceVertices.end());
    allSurfaceNormals.insert(allSurfaceNormals.end(), surfaceNormals.begin(),
                             surfaceNormals.end());
    auto nowSurfaceIndices = surfaceIndices;
    for (unsigned int& v : nowSurfaceIndices) v += i * surfaceNormals.size();

    allSurfaceIndices.insert(allSurfaceIndices.end(), nowSurfaceIndices.begin(),
                             nowSurfaceIndices.end());
  }

  // Generate a vertex array (VAO) and two vertex buffer objects (VBO).
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO_positions);
  glGenBuffers(1, &VBO_normals);

  // Bind to the VAO.
  glBindVertexArray(VAO);

  for (int i = 0; i < (int)positions.size(); ++i) {
    for (int j = 0; j < (int)surfaceVertices.size(); ++j) {
      allSurfaceVertices[i * surfaceVertices.size() + j] += positions[i];
    }
  }

  // Bind to the first VBO - We will use it to store the vertices
  glBindBuffer(GL_ARRAY_BUFFER, VBO_positions);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * allSurfaceVertices.size(),
               allSurfaceVertices.data(), GL_STATIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

  for (int i = 0; i < (int)positions.size(); ++i) {
    for (int j = 0; j < (int)surfaceVertices.size(); ++j) {
      allSurfaceVertices[i * surfaceVertices.size() + j] -= positions[i];
    }
  }

  // Bind to the second VBO - We will use it to store the normals
  glBindBuffer(GL_ARRAY_BUFFER, VBO_normals);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * allSurfaceNormals.size(),
               allSurfaceNormals.data(), GL_STATIC_DRAW);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

  // Generate EBO, bind the EBO to the bound VAO and send the data
  glGenBuffers(1, &EBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               sizeof(unsigned int) * allSurfaceIndices.size(),
               allSurfaceIndices.data(), GL_STATIC_DRAW);

  // Unbind the VBOs.
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

Fluid::~Fluid() {
  // Delete the VBOs and the VAO.
  glDeleteBuffers(1, &VBO_positions);
  glDeleteBuffers(1, &EBO);
  glDeleteVertexArrays(1, &VAO);
}

void Fluid::Draw(const glm::mat4& viewProjMtx, GLuint shader) {
  // actiavte the shader program
  glUseProgram(shader);

  // get the locations and send the uniforms to the shader
  glUniformMatrix4fv(glGetUniformLocation(shader, "viewProj"), 1, false,
                     (float*)&viewProjMtx);
  glUniformMatrix4fv(glGetUniformLocation(shader, "model"), 1, GL_FALSE,
                     (float*)&kModel);
  glUniform3fv(glGetUniformLocation(shader, "DiffuseColor"), 1, &kColor[0]);

  // Bind the VAO
  glBindVertexArray(VAO);

  // draw the points using triangles, indexed with the EBO
  glDrawElements(GL_TRIANGLES, allSurfaceIndices.size(), GL_UNSIGNED_INT, 0);

  // Unbind the VAO and shader program
  glBindVertexArray(0);
  glUseProgram(0);
}

void Fluid::Update() {
  // Bind to the VAO.
  glBindVertexArray(VAO);

  for (int i = 0; i < (int)positions.size(); ++i) {
    for (int j = 0; j < (int)surfaceVertices.size(); ++j) {
      allSurfaceVertices[i * surfaceVertices.size() + j] += positions[i];
    }
  }

  // Bind to the first VBO - We will use it to store the vertices
  glBindBuffer(GL_ARRAY_BUFFER, VBO_positions);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * allSurfaceVertices.size(),
               allSurfaceVertices.data(), GL_STATIC_DRAW);

  for (int i = 0; i < (int)positions.size(); ++i) {
    for (int j = 0; j < (int)surfaceVertices.size(); ++j) {
      allSurfaceVertices[i * surfaceVertices.size() + j] -= positions[i];
    }
  }
}

glm::vec3 Fluid::uniformSampleInsideSphere(
    float radius, std::default_random_engine& generator) {
  std::normal_distribution<float> normal(0.0, 1.0);
  float x = normal(generator);
  float y = normal(generator);
  float z = normal(generator);
  float normalize = pow(x * x + y * y + z * z, 0.5);

  x /= normalize;
  y /= normalize;
  z /= normalize;

  std::uniform_real_distribution<float> uniform(0.0, 1.0);
  float k = radius * pow(uniform(generator), 1.f / 3.f);
  return k * glm::vec3(x, y, z);
}

void Fluid::CreateSphereSurface() {
  // Copy and modify from http://www.songho.ca/opengl/gl_sphere.html
  float x, y, z, xy;                                   // vertex position
  float nx, ny, nz, lengthInv = 1.0f / kRenderRadius;  // vertex normal

  float sectorStep = 2 * PI / (float)sectorCount;
  float stackStep = PI / (float)stackCount;
  float sectorAngle, stackAngle;

  for (int i = 0; i <= stackCount; ++i) {
    stackAngle = PI / 2 - i * stackStep;    // starting from pi/2 to -pi/2
    xy = kRenderRadius * cosf(stackAngle);  // r * cos(u)
    z = kRenderRadius * sinf(stackAngle);   // r * sin(u)

    // add (sectorCount+1) vertices per stack
    // the first and last vertices have same position and normal, but different
    // tex coords
    for (int j = 0; j <= sectorCount; ++j) {
      sectorAngle = j * sectorStep;  // starting from 0 to 2pi

      // vertex position (x, y, z)
      x = xy * cosf(sectorAngle);  // r * cos(u) * cos(v)
      y = xy * sinf(sectorAngle);  // r * cos(u) * sin(v)
      surfaceVertices.emplace_back(x, y, z);

      // normalized vertex normal (nx, ny, nz)
      nx = x * lengthInv;
      ny = y * lengthInv;
      nz = z * lengthInv;
      surfaceNormals.emplace_back(nx, ny, nz);
    }
  }

  int k1, k2;
  for (int i = 0; i < stackCount; ++i) {
    k1 = i * (sectorCount + 1);  // beginning of current stack
    k2 = k1 + sectorCount + 1;   // beginning of next stack

    for (int j = 0; j < sectorCount; ++j, ++k1, ++k2) {
      // 2 triangles per sector excluding first and last stacks
      // k1 => k2 => k1+1
      if (i != 0) {
        surfaceIndices.push_back(k1);
        surfaceIndices.push_back(k2);
        surfaceIndices.push_back(k1 + 1);
      }

      // k1+1 => k2 => k2+1
      if (i != (stackCount - 1)) {
        surfaceIndices.push_back(k1 + 1);
        surfaceIndices.push_back(k2);
        surfaceIndices.push_back(k2 + 1);
      }
    }
  }
}

glm::vec3 Fluid::GetPosition(int i) { return positions[i]; }

int Fluid::GetNumOfParticles() { return positions.size(); }

void Fluid::SetPosition(int i, glm::vec3 p) { positions[i] = p; }

float Fluid::GetRadius() { return radius; }
