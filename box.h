#ifndef _BOX_H_
#define _BOX_H_

#include "core.h"
#include "mesh.h"

class Box : public Mesh {
 private:
  const glm::mat4 kModel = glm::mat4(1.0f);
  const glm::vec3 kColor = glm::vec3(1.0f);

  GLuint VAO;
  GLuint VBO_positions, EBO;

  std::vector<glm::vec3> positions;
  std::vector<unsigned int> indices;

 public:
  Box(glm::vec3 boxMin, glm::vec3 boxMax);

  ~Box();

  void Draw(const glm::mat4& viewProjMtx, GLuint shader);

  void Update();
};

#endif
