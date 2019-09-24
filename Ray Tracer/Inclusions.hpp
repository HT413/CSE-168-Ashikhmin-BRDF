// Global inclusions file
#ifndef GLOBAL_DEFS_H
#define GLOBAL_DEFS_H

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <vector>
#include <iostream>

// Global definitions
using namespace std;
using namespace glm;
#define PI 3.1415926536f
#define DEFAULT_LARGE_NUM 1000000000
#define MAX_DEPTH 5

// Global vars
const mat4 IDENTITY_MAT = mat4(1.f);


#endif