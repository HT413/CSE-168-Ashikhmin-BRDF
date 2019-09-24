// Definitions for objects used in the scene
#ifndef SCENE_INFO_HPP
#define SCENE_INFO_HPP

#include "SceneObject.hpp"
#include "Bitmap.h"

struct Light{
	float intensity = 1.f;
	vec4 position;
	vec3 attenuation = vec3(1, 0, 1);
	// Actual color: intensity * color
	vec3 color = vec3(1.f, 1.f, 1.f);
};

struct Scene{
	vector<Light> lights;
	vec3 backgroundColor;
	vector<SceneObject*> objects;

	void addObject(SceneObject* obj){ 
		objects.push_back(obj); 
	}
	void addLight(Light l){
		lights.push_back(l);
	}
	bool hitsObject(Ray&, vec4&);
};

extern Scene* currentScene;

struct Camera{
	vec3 cam_pos, cam_look_at, cam_up;
	int width, height;
	float fovy, fovx, aspect;
	mat4 view;
	Bitmap* output;
	int xSamples = 1, ySamples = 1;
	bool jittering = false;
	bool shirley = false;

	void setSuperSample(int x, int y){ xSamples = x; ySamples = y; }
	void setJitter(bool b){ jittering = b; }
	void setShirley(bool b){ shirley = b; }
	void resolution(int w, int h){
		width = w; height = h;
		aspect = float(w) / float(h);
		output = new Bitmap(w, h);
	}
	void fov(float fov){
		fovy = fov * PI / 180.f;
		fovx = 2 * atanf(tanf(fovy / 2.f) * width / height);
	}
	void lookAt(){
		vec3 d = cam_pos;
		vec3 c = normalize(d - cam_look_at);
		vec3 a = normalize(cross(cam_up, c));
		vec3 b = cross(c, a);
		view = mat4(a[0], a[1], a[2], 0, b[0], b[1], b[2], 0, c[0], c[1], c[2], 0, d[0], d[1], d[2], 1);
	}
	void lookAt(vec3& e, vec3& d, vec3& up){
		cam_pos = e; cam_look_at = d; cam_up = up;
		lookAt();
	};
	void saveBitmap(const char* filename){
		if(output->SaveBMP(filename)){
			cout << "Successfully saved BMP file with name " << filename << " of dimensions "
				 << width << "x" << height << "!" << endl;
		}
		else{
			cerr << "Could not save a BMP file with name " << filename << "...\n";
		}
		delete(output);
	}
	void render(Scene);
	vec3 calcColor(Scene, float, float);
};

#endif