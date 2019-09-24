// Main file for CSE 168 projects

#include "Inclusions.hpp"
#include "Scene.hpp"
#include <chrono>

void runProject1();
void runProject2();
void runProject3();
void runProject4();
Scene* currentScene; // The scene we're working with

// Main func. Simply runs project 1 and creates a bitmap out of it.
void main(){
	using namespace chrono;
	milliseconds startTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
	cout << "Running project. Please be patient." << endl;
	runProject4();
	milliseconds endTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
	cout << "Total execution time: " << float((endTime - startTime).count()) / 1000.f << " seconds." << endl;
	cout << "Finished running project. Press Enter to quit." << endl;
	getchar();
}

// Project 1 scene generation
void runProject1(){
	// Create the scene
	Scene scene;
	currentScene = &scene;
	scene.backgroundColor = vec3(.8f, .9f, 1.f);

	// Now create the boxes
	Box* box1 = new Box(5.f, .1f, 5.f);
	scene.addObject(box1);
	
	// Instance 1 of box 2
	Box* box2_1 = new Box(1.f, 1.f, 1.f);
	box2_1->setModel(rotate(IDENTITY_MAT, 0.5f, vec3(1, 0, 0)));
	box2_1->set(3, 1, 1.f);
	scene.addObject(box2_1);
	
	// Instance 2 of box 2
	Box* box2_2 = new Box(1.f, 1.f, 1.f);
	box2_2->setModel(rotate(IDENTITY_MAT, 1.f, vec3(0, 1, 0)));
	box2_2->set(3, vec4(-1, 0, 1, 1));
	scene.addObject(box2_2);
	
	/**/// Now create the lights
	Light sunLight;
	sunLight.color = vec3(1.f, 1.f, .9f);
	sunLight.intensity = .5f;
	sunLight.position = vec4(-.5f, -1.f, -.5f, 0.f);
	scene.addLight(sunLight);
	/**/
	Light redLight;
	redLight.color = vec3(1.f, .2f, .2f);
	redLight.intensity = 2.f;
	redLight.position = vec4(2.f, 2.f, 0.f, 1.f);
	scene.addLight(redLight);
	/**/
	// Create the camera
	Camera cam;
	cam.lookAt(vec3(2, 2, 5), vec3(0, 0, 0), vec3(0, 1, 0));
	cam.resolution(800, 600);
	cam.aspect = 1.33f;
	cam.fov(40.0f);

	// Render image
	cam.render(scene);
	cam.saveBitmap("project1.bmp");
}

// Project 2 scene generation
void runProject2() {
	// Create scene
	Scene scn;
	currentScene = &scn;
	scn.backgroundColor = vec3(0.8f, 0.8f, 1.0f);

	// Create ground
	SceneObject * ground = new Box(5.f, .1f, 5.f);
	scn.addObject(ground);

	// Create dragon
	PLY * dragon = new PLY("dragon.ply");
	dragon->smooth();

	BoxTreeObject * tree = new BoxTreeObject();
	tree->Construct(dragon);
	scn.addObject(tree);

	// Create instance
	BoxTreeObject * inst = new BoxTreeObject();
	inst->Construct(dragon);
	//InstanceObject * inst = new InstanceObject(tree);
	glm::mat4x4 mtx = eulerAngleY(PI);
	mtx[3] = glm::vec4(-0.05f, 0.0f, -0.1f, 1.0f);
	inst->setModel(mtx);
	//scn.addObject(inst);

	// Create lights
	Light sunlgt;
	sunlgt.color = vec3(1.0f, 1.0f, 0.9f);
	sunlgt.intensity = 1.0f;
	sunlgt.position = vec4(2.0f, -3.0f, -2.0f, 0.f);
	scn.addLight(sunlgt);

	Light redlgt;
	redlgt.color = vec3(1.0f, 0.2f, 0.2f);
	redlgt.intensity = 0.02f;
	redlgt.position = vec4(-0.2f, 0.2f, 0.2f, 1.f);
	scn.addLight(redlgt);

	Light bluelgt;
	bluelgt.color = vec3(0.2f, 0.2f, 1.0f);
	bluelgt.intensity = 0.02f;
	bluelgt.position = vec4(0.1f, 0.1f, 0.3f, 1.f);
	scn.addLight(bluelgt);

	// Create camera
	Camera cam;
	cam.lookAt(vec3(-0.1f, 0.1f, 0.2f), vec3(-0.05f, 0.12f, 0.0f), vec3(0, 1.0f, 0));
	cam.resolution(800, 600);
	cam.aspect = 1.33f;
	cam.fov(40.0f);
	
	// Render image
	cam.render(scn);
	cam.saveBitmap("project2.bmp");
}

// Project 3 scene generation
void runProject3() {
	// Create scene
	Scene scn;
	currentScene = &scn;
	scn.backgroundColor = vec3(.8f, .9f, 1.f);

	// Create ground
	LambertMaterial * groundMtl = new LambertMaterial(0.25f, 0.25f, 0.25f);

	SceneObject * ground = new Box(2.f, .11f, 2.f);
	ground->setMaterial(groundMtl);
	scn.addObject(ground);

	// Load dragon mesh
	PLY * dragon = new PLY("dragon.ply");

	// Create box tree
	//BoxTreeObject * tree = new BoxTreeObject();
	//tree->Construct(dragon);

	// Materials
	LambertMaterial * white = new LambertMaterial(0.7f, 0.7f, 0.7f);
	LambertMaterial * red = new LambertMaterial(0.7f, 0.1f, 0.1f);
	MetalMaterial * metal = new MetalMaterial(0.95f, 0.64f, 0.54f);

	const int numDragons=4;
	Material * mtl[numDragons]={white, metal, red, white};
	// Create dragon instances
	glm::mat4 mtx;
	for(int i=0; i<numDragons; i++) {
		mtx[3]=glm::vec4(0.0f, 0.0f, 0.3f*(float(i)/float(numDragons-1)-0.5f), 1.0f);
		BoxTreeObject * inst = new BoxTreeObject();
		inst->Construct(dragon);
		inst->setModel(mtx);
		inst->setMaterial(mtl[i]);
		scn.addObject(inst);
	}
	// Create lights
	Light sunlgt;
	sunlgt.color = vec3(1.0f, 1.0f, 0.9f);
	sunlgt.intensity = 1.0f;
	sunlgt.position = vec4(2.0f, -3.0f, -2.0f, 0.f);
	scn.addLight(sunlgt);

	// Create camera
	Camera cam;
	cam.lookAt(vec3(-0.5f, 0.25f, -0.2f), vec3(0.0f, 0.15f, 0.0f), vec3(0, 1.0f, 0));
	cam.resolution(640, 480);
	cam.aspect = 1.33f;
	cam.fov(40.0f);
	cam.setSuperSample(2, 2);
	cam.setJitter(true);
	cam.setShirley(true);

	// Render image
	cam.render(scn);
	cam.saveBitmap("project3.bmp");
}

// Project 4 scene generation
void runProject4() {
	// Create scene
	Scene scn;
	currentScene = &scn;
	scn.backgroundColor = vec3(.8f, .9f, 1.f);
	// Materials
	const int nummtls=4;
	AshikhminMaterial * mtl[nummtls];
	// Diffuse
	mtl[0] = new AshikhminMaterial();
	mtl[0]->setSpecularLevel(0.0f);
	mtl[0]->setDiffuseLevel(1.0f);
	mtl[0]->setDiffuseColor(vec3(0.7f, 0.7f, 0.7f));
	// Roughened copper
	mtl[1] = new AshikhminMaterial();
	mtl[1]->setDiffuseLevel(0.0f);
	mtl[1]->setSpecularLevel(1.0f);
	mtl[1]->setSpecularColor(vec3(0.9f, 0.6f, 0.5f));
	mtl[1]->setRoughness(100.0f, 100.0f);
	// Anisotropic gold
	mtl[2] = new AshikhminMaterial();
	mtl[2]->setDiffuseLevel(0.0f);
	mtl[2]->setSpecularLevel(1.0f);
	mtl[2]->setSpecularColor(vec3(0.95f, 0.7f, 0.3f));
	mtl[2]->setRoughness(1.0f, 1000.0f);
	// Red plastic
	mtl[3] = new AshikhminMaterial();
	mtl[3]->setDiffuseColor(vec3(1.0f, 0.1f, 0.1f));
	mtl[3]->setDiffuseLevel(0.8f);
	mtl[3]->setSpecularLevel(0.2f);
	mtl[3]->setSpecularColor(vec3(1.0f, 1.0f, 1.0f));
	mtl[3]->setRoughness(1000.0f, 1000.0f);
	// Load dragon mesh
	PLY * dragon = new PLY("dragon.ply");
	// Create box tree
	BoxTreeObject tree;
	tree.Construct(dragon);
	// Create dragon instances
	glm::mat4 mtx;
	for(int i=0; i<nummtls; i++) {
		mtx[3] = vec4(0.0f, 0.0f, -0.1f*float(i), 1.0f);
		BoxTreeObject * inst = new BoxTreeObject();
		inst->Construct(dragon);
		inst->setModel(mtx);
		inst->setMaterial(mtl[i]);
		scn.addObject(inst);
	}	
	// Create ground
	LambertMaterial * groundMtl = new LambertMaterial(0.3f, 0.3f, 0.35f);

	SceneObject * ground = new Box(2.f, .11f, 2.f);
	ground->setMaterial(groundMtl);
	scn.addObject(ground);
	// Create lights
	Light sunlgt;
	sunlgt.color = vec3(1.0f, 1.0f, 0.9f);
	sunlgt.intensity = 1.0f;
	sunlgt.position = vec4(2.0f, -3.0f, -2.0f, 0.f);
	scn.addLight(sunlgt);	
	// Create camera
	Camera cam;
	cam.lookAt(vec3(-0.5f, 0.25f, -0.2f), vec3(0.0f, 0.15f, -0.15f), vec3(0, 1.0f, 0));
	cam.resolution(800, 600);
	cam.aspect = 1.33f;
	cam.fov(40.0f);
	cam.setSuperSample(2, 2);
	cam.setJitter(true);
	cam.setShirley(true);

	// Render image
	cam.render(scn);
	cam.saveBitmap("project4.bmp");
}