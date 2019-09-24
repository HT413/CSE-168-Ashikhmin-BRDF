// Definitions for scene objects
#ifndef SCENE_OBJECT_HPP
#define SCENE_OBJECT_HPP

#include "Inclusions.hpp"

struct Ray{
	vec3 origin, direction;
};

struct Material{
protected:
	vec3 color = vec3(0, 0, 0);
public:
	vec3& getColor(){ return color; }
	virtual vec3 computeColor(vec3&, vec3&, vec3&, vec3&) = 0;
	virtual void generateSample(const vec3&, const vec3&, const vec3&, vec3&, vec3&) = 0;
	virtual void generateSample(const vec3&, const vec3&, const vec3&, const vec3&, vec3&, vec3&) = 0;
};

struct LambertMaterial : public Material{
	LambertMaterial(float r, float g, float b){ color = vec3(r, g, b); }
	void generateSample(const vec3&, const vec3&, const vec3&, vec3&, vec3&);
	void generateSample(const vec3& a, const vec3& b, const vec3& c, const vec3& d, vec3& e, vec3& f){
		generateSample(a, b, c, e, f);
	}
	vec3 computeColor(vec3& normal, vec3& dir, vec3& inColor, vec3& lDir);
};

struct MetalMaterial: public Material{
	MetalMaterial(float r, float g, float b){ color = vec3(r, g, b); }
	void generateSample(const vec3&, const vec3&, const vec3&, vec3&, vec3&);
	void generateSample(const vec3& a, const vec3& b, const vec3& c, const vec3& d, vec3& e, vec3& f){
		generateSample(a, b, c, e, f);
	}
	vec3 computeColor(vec3& normal, vec3& dir, vec3& inColor, vec3& lDir);
};

struct AshikhminMaterial: public Material{
private:
	float diffuseLevel = 0, specularLevel = 0;
	float n_u = 0, n_v = 0;
protected:
	vec3 specularColor = vec3(0, 0, 0);
public:
	void setDiffuseLevel(float f){ diffuseLevel = f; };
	void setSpecularLevel(float f){ specularLevel = f; };
	void setDiffuseColor(vec3 v){ color = v; };
	void setSpecularColor(vec3 v){ specularColor = v; }
	void setRoughness(float f1, float f2){ n_u = f1; n_v = f2; };
	void generateSample(const vec3&, const vec3&, const vec3&, vec3&, vec3&);
	vec3 computeColor(vec3& normal, vec3& dir, vec3& inColor, vec3& lDir);
	void generateSample(const vec3&, const vec3&, const vec3&, const vec3&, vec3&, vec3&) override;
};

vec3 intersectTriangle(Material*, vec3&, Ray&, vec3&, vec3&, vec3&, float&, int);
vec3 intersectTriangle(Material*, vec3&, Ray&, vec3&, vec3&, vec3&, float&, vec3&, int);
bool intersectsTriangle(Ray&, vec3&, vec3&, vec3&, vec4&);

class SceneObject{
protected:
	vector<vec3> positions;
	vector<vec3> normals;
	vector<int> indices;
	Material* mtl;
	mat4 model = mat4(1.f);
	mat4 invModel = mat4(1.f);
	mat4 parentModel = mat4(1.f);

public:
	vector<vec3>& getPositions(){ return positions; }
	vector<int>& getIndices(){ return indices; }
	mat4 compositeModel(){ return parentModel * model; }
	vec3& vertexAt(int i){ return positions[i]; }
	int& indexAt(int i){ return indices[i]; }
	vec3& normalAt(int i){ return normals[i]; }
	mat4 getParentModel(){ return parentModel; }
	void setParentModel(mat4 m){ parentModel = m; }
	void setModel(mat4&);
	void set(int, vec4&);
	void set(int, int, float);
	virtual vec4 intersection(vec3, Ray, float&, int) = 0;
	Material* getMaterial(){ return mtl; }
	void setMaterial(Material * m){ mtl = m; }
	virtual bool intersects(Ray&, vec4&) = 0;
};

class InstanceObject: public SceneObject{
private:
	SceneObject* obj;
public:
	InstanceObject(SceneObject* o){ obj = o; }
	vec4 intersection(vec3 v, Ray r, float& f, int d) override{
		mat4& prev = obj->getParentModel();
		Material* prevMat = obj->getMaterial();
		obj->setParentModel(model);
		obj->setMaterial(mtl);
		vec4 result = obj->intersection(v, r, f, d); 
		obj->setParentModel(prev);
		obj->setMaterial(prevMat);
		return result;
	}
	bool intersects(Ray& r, vec4& max) override{
		mat4& prev = obj->getParentModel();
		Material* prevMat = obj->getMaterial();
		obj->setParentModel(model);
		obj->setMaterial(mtl);
		bool result = obj->intersects(r, max);
		obj->setParentModel(prev);
		obj->setMaterial(prevMat);
		return result;
	};
};

class Box:public SceneObject{
public:
	Box(float width, float height, float depth){
		//type = box;
		// Generate the positions for this box
		float w = width * .5f;
		float h = height * .5f;
		float d = depth * .5f;

		// Push the corners
		// Front face, start bottom left, go ccw
		positions.push_back(vec3(-w, -h, d));
		positions.push_back(vec3(w, -h, d));
		positions.push_back(vec3(w, h, d));
		positions.push_back(vec3(-w, h, d));
		// Back face, start bottom left, go ccw
		positions.push_back(vec3(-w, -h, -d));
		positions.push_back(vec3(w, -h, -d));
		positions.push_back(vec3(w, h, -d));
		positions.push_back(vec3(-w, h, -d));

		// Now push the indices
		indices.assign({0, 1, 2, 0, 2, 3, // Front face
		                4, 6, 5, 4, 7, 6, // Back face
		                0, 5, 1, 0, 4, 5, // Bottom face
		                3, 2, 6, 3, 6, 7, // Top face
		                1, 5, 6, 1, 6, 2, // Right face
		                0, 7, 4, 0, 3, 7}); // Left face
	}
	vec4 intersection(vec3, Ray, float&, int) override;
	bool intersects(Ray&, vec4&) override;
};

class PLY: public SceneObject{
private:
	int numVertices;

public:
	vector<int> triangles;
	float minX = INFINITY, minY = INFINITY, minZ = INFINITY;
	float maxX = -INFINITY, maxY = -INFINITY, maxZ = -INFINITY;

	PLY(const char* filename){
		loadPLY(filename);
	};

	bool loadPLY(const char*);
	void smooth();
	vec4 intersection(vec3, Ray, float&, int) { return vec4(0, 0, 0, 0); };
	bool intersects(Ray&, vec4&) override { return false; }
};

class BoxTreeObject;

struct TreeNode{
	PLY* o;
	BoxTreeObject* tree;
	TreeNode *left = nullptr, *right = nullptr;
	TreeNode(PLY* p, BoxTreeObject* bt){ 
		o = p; tree = bt; 
		min = vec3(DEFAULT_LARGE_NUM, DEFAULT_LARGE_NUM, DEFAULT_LARGE_NUM);
		max = vec3(-DEFAULT_LARGE_NUM, -DEFAULT_LARGE_NUM, -DEFAULT_LARGE_NUM);
	}
	vec3 min, max;
	vector<int> triangleIndices;
	void Construct(PLY*, BoxTreeObject*, vector<int> indices);
	vec4 intersection(vec3&, Ray, float&, int);
	float intersect(Ray);
	bool intersects(Ray&, vec4&);
};

class BoxTreeObject: public SceneObject{
private:
	TreeNode *root;

public:
	BoxTreeObject(){}
	void Construct(PLY*);
	vec4 intersection(vec3, Ray, float&, int) override;
	bool intersects(Ray&, vec4&) override;
};

#endif