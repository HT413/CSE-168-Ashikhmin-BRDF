#include "Scene.hpp"
#include <future>

int colorInt(vec4 c){
	int r = (c.x <= 0.f)? 0 : (c.x >= 1.f)? 255 : int(c.x * 256.f);
	int g = (c.y <= 0.f)? 0 : (c.y >= 1.f)? 255 : int(c.y * 256.f);
	int b = (c.z <= 0.f)? 0 : (c.z >= 1.f)? 255 : int(c.z * 256.f);
	int a = int(c.a * 256.f);
	//return (b << 24) | (r << 16) | (g << 8) | a;
	return (r << 16) | (g << 8) | b;
}

float randomF(float maxAbs){
	return ((rand() % 32768) * maxAbs / 16383.5f) - maxAbs;
}

// Render func. Loops through every pixel to render
void Camera::render(Scene s){
	size_t max = width * height;
	size_t cores = thread::hardware_concurrency() - 1;
	volatile atomic<size_t> count(0);
	vector<future<void>> future_vector;
	cout << "Beginning ray trace. Please be patient..." << endl;
	while(cores--){
		future_vector.emplace_back(
			async([=, &count](){
			while(true){
				size_t index = count++;
				if(index >= max) break;
				int x = index % width;
				int y = index / width;
				if(((y * 10) % height == 0) && (x == 0))
					cout << (y * 100 / height) << "% done!" << endl;
				vec3 sum(0, 0, 0);
				for(int subX = 0; subX < xSamples; subX++){
					for(int subY = 0; subY < ySamples; subY++){
						float subCenterX = (subX + 0.5f) / float(xSamples);
						float subCenterY = (subY + 0.5f) / float(ySamples);
						if(jittering){
							subCenterX += randomF(0.5f / xSamples);
							subCenterY += randomF(0.5f / xSamples);
						}
						if(shirley){
							subCenterX = (subCenterX < .5f)? -.5f + sqrtf(2 * subCenterX) :
								1.5f - sqrtf(2 - 2 * subCenterX);
							subCenterY = (subCenterY < .5f)? -.5f + sqrtf(2 * subCenterY) :
								1.5f - sqrtf(2 - 2 * subCenterY);
						}
						sum += calcColor(s, float(x) + subCenterX, float(y) + subCenterY);
					}
				}
				sum = sum * (1.f / xSamples / ySamples);
				output->SetPixel(x, y, colorInt(vec4(sum, 1)));
			}
		}));
	}
}

vec3 Camera::calcColor(Scene s, float x, float y){
	// Generate ray to this particular pixel
	Ray r;
	r.origin = cam_pos;
	float alpha = tan(fovx / 2.f) * (x - float(width / 2)) / float(width / 2);
	float beta = tan(fovy / 2.f) * (y - float(height / 2)) / float(height / 2);
	r.direction = normalize(alpha * vec3(view[0]) + beta * vec3(view[1]) - vec3(view[2]));
	// Now loop through the objects/triangles and test
	float t = DEFAULT_LARGE_NUM; // Default: very large timestep
	vec4 finalCol = vec4(s.backgroundColor, 1.f);
	for(SceneObject *obj: s.objects){
		finalCol = obj->intersection(vec3(finalCol), r, t, 0);
	}
	if(finalCol.x > 1) finalCol.x = 1;
	if(finalCol.y > 1) finalCol.y = 1;
	if(finalCol.z > 1) finalCol.z = 1;
	return finalCol;
}

bool Scene::hitsObject(Ray& r, vec4& lPos){
	bool result = false;
	for(SceneObject *o: objects){
		result = o->intersects(r, lPos);
		if(result) break;
	}
	return result;
}