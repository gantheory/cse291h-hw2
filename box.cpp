#include "box.h"
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <random>
#include "core.h"

#define PI 3.1415926f

Box::Box(glm::vec3 boxMin, glm::vec3 boxMax) : boxMin(boxMin), boxMax(boxMax) {
  std::vector<glm::vec3> points = {boxMin, boxMax};
  positions = {// [0, 0, 0] <-> [0, 0, 1]
               glm::vec3(points[0].x, points[0].y, points[0].z),
               glm::vec3(points[0].x, points[0].y, points[1].z),
               // [0, 0, 0] <-> [1, 0, 0]
               glm::vec3(points[0].x, points[0].y, points[0].z),
               glm::vec3(points[1].x, points[0].y, points[0].z),
               // [1, 0, 0] <-> [1, 0, 1]
               glm::vec3(points[1].x, points[0].y, points[0].z),
               glm::vec3(points[1].x, points[0].y, points[1].z),
               // [0, 0, 1] <-> [1, 0, 1]
               glm::vec3(points[0].x, points[0].y, points[1].z),
               glm::vec3(points[1].x, points[0].y, points[1].z),
               // [0, 0, 1] <-> [0, 1, 1]
               glm::vec3(points[0].x, points[0].y, points[1].z),
               glm::vec3(points[0].x, points[1].y, points[1].z),
               // [0, 0, 0] <-> [0, 1, 0]
               glm::vec3(points[0].x, points[0].y, points[0].z),
               glm::vec3(points[0].x, points[1].y, points[0].z),
               // [1, 0, 0] <-> [1, 1, 0]
               glm::vec3(points[1].x, points[0].y, points[0].z),
               glm::vec3(points[1].x, points[1].y, points[0].z),
               // [1, 0, 1] <-> [1, 1, 1]
               glm::vec3(points[1].x, points[0].y, points[1].z),
               glm::vec3(points[1].x, points[1].y, points[1].z),
               // [0, 1, 0] <-> [0, 1, 1]
               glm::vec3(points[0].x, points[1].y, points[0].z),
               glm::vec3(points[0].x, points[1].y, points[1].z),
               // [0, 1, 0] <-> [1, 1, 0]
               glm::vec3(points[0].x, points[1].y, points[0].z),
               glm::vec3(points[1].x, points[1].y, points[0].z),
               // [1, 1, 0] <-> [1, 1, 1]
               glm::vec3(points[1].x, points[1].y, points[0].z),
               glm::vec3(points[1].x, points[1].y, points[1].z),
               // [0, 1, 1] <-> [1, 1, 1]
               glm::vec3(points[0].x, points[1].y, points[1].z),
               glm::vec3(points[1].x, points[1].y, points[1].z)};
  indices.resize(positions.size());
  iota(indices.begin(), indices.end(), 0);

  // Generate a vertex array (VAO) and two vertex buffer objects (VBO).
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO_positions);

  // Bind to the VAO.
  glBindVertexArray(VAO);

  // Bind to the first VBO - We will use it to store the vertices
  glBindBuffer(GL_ARRAY_BUFFER, VBO_positions);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * positions.size(),
               positions.data(), GL_STATIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

  // Generate EBO, bind the EBO to the bound VAO and send the data
  glGenBuffers(1, &EBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices.size(),
               indices.data(), GL_STATIC_DRAW);

  // Unbind the VBOs.
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

Box::~Box() {
  // Delete the VBOs and the VAO.
  glDeleteBuffers(1, &VBO_positions);
  glDeleteBuffers(1, &EBO);
  glDeleteVertexArrays(1, &VAO);
}

void Box::Draw(const glm::mat4& viewProjMtx, GLuint shader) {
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
  glDrawElements(GL_LINES, indices.size(), GL_UNSIGNED_INT, 0);

  // Unbind the VAO and shader program
  glBindVertexArray(0);
  glUseProgram(0);
}

void Box::Update() {}

glm::vec3 Box::GetBoxMin() { return boxMin; }

glm::vec3 Box::GetBoxMax() { return boxMax; }
