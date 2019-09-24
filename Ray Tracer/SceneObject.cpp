#include "SceneObject.hpp"
#include "Scene.hpp"

const int MAX_PER_BOX = 20;

float randomPF(float max){
	return (rand() % 32768) / 32767.f * max;
}

void SceneObject::setModel(mat4& m){
	model = m;
	invModel = inverse(model);
}

void SceneObject::set(int col, int row, float val){
	model[col][row] = val;
	invModel = inverse(model);
}

void SceneObject::set(int col, vec4& v){
	model[col] = v;
	invModel = inverse(model);
}

vec3 LambertMaterial::computeColor(vec3& normal, vec3& dir, vec3& inColor, vec3& lDir){
	if(dot(lDir, normal) < 0.f) return vec3(0, 0, 0);
	return color * inColor;
}

vec3 MetalMaterial::computeColor(vec3& normal, vec3& dir, vec3& inColor, vec3& lDir){
	return vec3(0, 0, 0);
}

void LambertMaterial::generateSample(const vec3& intersection, const vec3& norm, const vec3& in, vec3& out, vec3& col){
	vec3 rand = normalize(vec3(randomPF(1.f), randomPF(1.f), randomPF(1.f)));
	vec3 ivec = normalize(cross(norm, rand));
	vec3 kvec = normalize(cross(ivec, norm));
	float u = 2 * PI * randomPF(1.f);
	float t = randomPF(1.f);
	float v = sqrtf(1 - t);
	out = normalize(v * cosf(u) * ivec + sqrtf(t) * norm + v * sinf(u) * kvec);
	col = color;
}

void MetalMaterial::generateSample(const vec3& intersection, const vec3& norm, const vec3& in, vec3& out, vec3& col){
	out = (dot(in, norm) < 0.f)? in + 2 * dot(-in, norm) * norm : -in + 2 * dot(in, norm) * norm;
	out = normalize(out);
	// Assume total reflection for now
	col = color;
}

vec3 AshikhminMaterial::computeColor(vec3& normal, vec3& dir, vec3& inColor, vec3& lDir){
	vec3 halfVec = normalize(dir + lDir);
	float Fresnel = specularLevel + (1.f - specularLevel) * powf((1.f - dot(dir, halfVec)), 5);
	
	vec3 uVec = cross(vec3(0, 1, 0), normal);
	if(length(uVec) < .0001f) uVec = cross(vec3(1, 0, 0), normal);
	vec3 vVec = normalize(cross(normal, uVec));
	
	float specular = sqrtf((n_u + 1) * (n_v + 1)) / (8.f /* PI */) * Fresnel 
		* powf(dot(normal, halfVec), 
		((n_u * powf(dot(halfVec, uVec), 2) + (n_v * powf(dot(halfVec, vVec), 2)))
			/ (1 - powf(dot(halfVec, normal), 2))))
		/ (dot(halfVec, dir) * max(dot(normal, dir), dot(normal, lDir)));
	float diffuse = 28.f * diffuseLevel / (23.f /* PI */) * (1.f - specularLevel) *
		(1.f - powf(1.f - dot(normal, dir) / 2.f, 5)) *
		(1.f - powf(1.f - dot(normal, lDir) / 2.f, 5));
	return inColor * (diffuse * color + specular * specularColor);
}

void AshikhminMaterial::generateSample(const vec3& intersection, const vec3& norm, const vec3& in, const vec3& lDir, vec3& out, vec3& col){
	float which = randomPF(1.f);
	if(which > specularLevel){
		vec3 rand = normalize(vec3(randomPF(1.f), randomPF(1.f), randomPF(1.f)));
		vec3 ivec = normalize(cross(norm, rand));
		vec3 kvec = normalize(cross(ivec, norm));
		float u = 2 * PI * randomPF(1.f);
		float t = randomPF(1.f);
		float v = sqrtf(1 - t);
		out = normalize(v * cosf(u) * ivec + sqrtf(t) * norm + v * sinf(u) * kvec);
		col = diffuseLevel * color;
	}
	else{
		float zeta1 = randomPF(1.f), zeta2 = randomPF(1.f);
		float phi = atanf(sqrtf((n_u + 1.f) / (n_v + 1.f)) * tanf(PI * zeta1 / 2.f));
		int quadrant = rand() % 4;
		if(quadrant == 1) phi = PI - phi;
		else if(quadrant == 2) phi = PI + phi;
		else if(quadrant == 3) phi = 2.f * PI - phi;

		float cosTheta = powf((1.f - zeta2), 1.f /
			(1.f + (n_u * powf(cosf(phi), 2)) + (n_v * powf(sinf(phi), 2))));
		float sinTheta = sqrtf(1.f - cosTheta * cosTheta);

		// Tangent vectors
		vec3 uVec = cross(vec3(0, 1, 0), norm);
		if(length(uVec) < .0001f) uVec = cross(vec3(1, 0, 0), norm);
		vec3 vVec = normalize(cross(norm, uVec));

		vec3 hVec = vec3(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
		hVec = normalize(hVec.y * uVec + hVec.z * norm + hVec.x * vVec);
		
		vec3 k2 = normalize(in + 2.f * dot(hVec, -in) * hVec);
		
		//if(dot(k2, norm) < 0) k2 = -k2;
		//vec3 k2 = (dot(in, hVec) < 0)? normalize(-in + 2.f * dot(hVec, in) * hVec)
		//		: normalize(in + 2.f * dot(hVec, -in) * hVec);

		//float p_h = sqrtf((n_u + 1.f) * (n_v + 1.f) / (2 * PI)) * powf(dot(norm, hVec),
		//	n_u * powf(cosf(phi), 2) + n_v * powf(sinf(phi), 2));
		//float pk2 = (dot(k2, norm) < 0)? 0.f : p_h / (4.f * dot(-in, hVec));
		
		//out = normalize(cosTheta * uVec + cosf(phi) * sinTheta * norm
		//	+ sinf(phi) * sinTheta * vVec);

		out = k2;
		col = specularLevel * specularColor;
	}
}

void AshikhminMaterial::generateSample(const vec3& intersection, const vec3& norm, const vec3& in, vec3& out, vec3& col){
	out = (dot(in, norm) < 0.f)? in + 2 * dot(-in, norm) * norm : -in + 2 * dot(in, norm) * norm;
	// Assume total reflection for now
	col = color;
}

vec3 intersectTriangle(Material* mtl, vec3& bg, Ray& r, vec3& a, vec3& b, vec3& c, float& t, vec3& nor, int depth){
	if(depth > MAX_DEPTH) return vec3(0, 0, 0);
	// Get the triangle normal
	vec3 normal = cross(b - a, c - a);
	vec3 n = normalize(normal);
	// Plane intersection
	if(dot(r.direction, n) == 0.f) return bg;
	float time = (dot((a - r.origin), n) / dot(r.direction, n));
	// Now check if intersection point is in triangle or not, if closer
	if(t < time || time <= 0.f) return bg;
	vec3 p = r.origin + time * r.direction;
	float det = dot(-r.direction, normal);
	if(det == 0.f) return bg;
	float alpha = dot(-r.direction, cross(p - a, c - a)) / det;
	float beta = dot(-r.direction, cross(b - a, p - a)) / det;
	if(!(alpha > 0.f && beta > 0.f && alpha + beta < 1.f)) return bg;
	t = time;
	// Calculate lighting and return that color
	vec3 finalCol = vec3(0, 0, 0);
	for(Light l: currentScene->lights){
		// Get the light direction
		vec3 direction;
		vec3 lightCol = l.intensity * l.color;
		// Directional light
		if(l.position.w < 0.001f){
			direction = normalize(vec3(-l.position));
		}
		// Point light
		else{
			direction = normalize(vec3(l.position) - p);
			lightCol *= (1.f / powf(length(vec3(l.position) - p), 2));
		}
		float nDotL = dot(nor, direction);
		if(depth == 0){
			// Compute diffuse lighting and add it
			float shadowFactor = 1.f;
			if(nDotL > 0.f){
				Ray rnew;
				rnew.direction = direction;
				rnew.origin = p + 0.001f * direction;
				if(currentScene->hitsObject(rnew, l.position))
					shadowFactor = 0.f;
			}
			finalCol += max(nDotL, 0.f) * mtl->computeColor(nor, -r.direction, lightCol, direction) * shadowFactor;
		}

		if(MAX_DEPTH > 0){
			Ray outgoing;
			vec3 outColor;
			float newTime = DEFAULT_LARGE_NUM;
			mtl->generateSample(p, nor, r.direction, direction, outgoing.direction, outColor);
			outgoing.origin = p + 0.005f * outgoing.direction;
			vec4 finalColor;
			vec3 background = vec3(currentScene->backgroundColor);
			for(SceneObject *obj: currentScene->objects){
				background = vec3(obj->intersection(background, outgoing, newTime, depth + 1));
			}
			finalColor = vec4(outColor * background, 1.f);
			finalCol += vec3(finalColor);
		}
	}
	return finalCol;
}

vec3 intersectTriangle(Material* mtl, vec3& bg, Ray& r, vec3& a, vec3& b, vec3& c, float& t, int depth){
	if(depth > MAX_DEPTH) return vec3(0, 0, 0);
	// Get the triangle normal
	vec3 normal = cross(b - a, c - a);
	vec3 n = normalize(normal);
	// Plane intersection
	if(dot(r.direction, n) == 0.f) return bg;
	float time = (dot((a - r.origin), n) / dot(r.direction, n));
	// Now check if intersection point is in triangle or not, if closer
	if(t < time || time <= 0.f) return bg;
	vec3 p = r.origin + time * r.direction;
	float det = dot(-r.direction, normal);
	if(det == 0.f) return bg;
	float alpha = dot(-r.direction, cross(p - a, c - a)) / det;
	float beta = dot(-r.direction, cross(b - a, p - a)) / det;
	if(!(alpha > 0.f && beta > 0.f && alpha + beta < 1.f)) return bg;
	t = time;
	// Calculate lighting and return that color
	vec3 finalCol = vec3(0, 0, 0);
	for(Light l: currentScene->lights){
		// Get the light direction
		vec3 direction;
		vec3 lightCol = l.intensity * l.color;
		// Directional light
		if(l.position.w < 0.001f){
			direction = normalize(vec3(-l.position));
		}
		// Point light
		else{
			direction = normalize(vec3(l.position) - p);
			lightCol *= (1.f / powf(length(vec3(l.position) - p), 2));
		}
		float nDotL = dot(n, direction);
		if(depth == 0){
			// Compute diffuse lighting and add it
			float shadowFactor = 1.f;
			if(nDotL > 0.f){
				Ray rnew;
				rnew.direction = direction;
				rnew.origin = p + 0.005f * direction;
				if(currentScene->hitsObject(rnew, l.position))
					shadowFactor = 0.f;
			}
			finalCol += max(nDotL, 0.f) * mtl->computeColor(n, -r.direction, lightCol, direction) * shadowFactor;
		}
	
		if(MAX_DEPTH > 0){
			Ray outgoing;
			vec3 outColor;
			float newTime = DEFAULT_LARGE_NUM;
			mtl->generateSample(p, n, r.direction, direction, outgoing.direction, outColor);
			outgoing.origin = p + 0.001f * outgoing.direction;
			vec4 finalColor;
			vec3 background = vec3(currentScene->backgroundColor);
			for(SceneObject *obj: currentScene->objects){
				background = vec3(obj->intersection(background, outgoing, newTime, depth + 1));
			}
			finalColor = vec4(outColor * background, 1.f);
			finalCol += vec3(finalColor);
		}
	}
	return finalCol;
}

bool intersectsTriangle(Ray& r, vec3& a, vec3& b, vec3& c, vec4& lPos){
	// Get the triangle normal
	vec3 normal = cross(b - a, c - a);
	vec3 n = normalize(normal);
	// Plane intersection
	if(dot(r.direction, n) == 0.f) return false;
	float time = (dot((a - r.origin), n) / dot(r.direction, n));
	// Now check if intersection point is in triangle or not, if closer
	if(time <= 0) return false;
	vec3 p = r.origin + time * r.direction;
	float det = dot(-r.direction, normal);
	if(det == 0.f) return false;
	float alpha = dot(-r.direction, cross(p - a, c - a)) / det;
	float beta = dot(-r.direction, cross(b - a, p - a)) / det;
	if(!(alpha > 0.f && beta > 0.f && alpha + beta < 1.f)) return false;
	if(lPos.w < 0.001f) return true;
	if(length(vec3(lPos) - r.origin) < length(p - r.origin)) return false;
	// Made it this far, it intersects!
	return true;
}

bool Box::intersects(Ray& r, vec4& max){
	mat4 compModel = compositeModel();
	bool result = false;
	for(unsigned int i = 0; i < indices.size(); i+=3){
		vec3 a = vec3(compModel * vec4(positions[indices[i]], 1.f));
		vec3 b = vec3(compModel * vec4(positions[indices[i+1]], 1.f));
		vec3 c = vec3(compModel * vec4(positions[indices[i+2]], 1.f));
		result = intersectsTriangle(r, a, b, c, max);
		if(result)
			break;
	}
	return result;
}

vec4 Box::intersection(vec3 bg, Ray r, float& t, int d){
	vec3 returnCol = bg;
	mat4 compModel = compositeModel();
	for(unsigned int i = 0; i < indices.size(); i+=3){
		vec3 a = vec3(compModel * vec4(positions[indices[i]], 1.f));
		vec3 b = vec3(compModel * vec4(positions[indices[i+1]], 1.f));
		vec3 c = vec3(compModel * vec4(positions[indices[i+2]], 1.f));
		returnCol = intersectTriangle(mtl, returnCol, r, a, b, c, t, d);
	}
	return vec4(returnCol, 1.f);
}

bool PLY::loadPLY(const char* filename){// Open file
	cout << "Loading file " << filename << endl;
	FILE *f=fopen(filename, "r");
	if(f==0) {
		printf("ERROR: MeshObject::LoadPLY()- Can't open '%s'\n", filename);
		return false;
	}
	// Read header
	char tmp[256];
	int numverts=0, numtris=0;
	int posprop=-99, normprop=-99;
	int props=0;
	while(1) {
		fgets(tmp, 256, f);
		if(strncmp(tmp, "element vertex", 14)==0)
			numverts=atoi(&tmp[14]);
		if(strncmp(tmp, "element face", 12)==0)
			numtris=atoi(&tmp[12]);
		if(strncmp(tmp, "property", 8)==0) {
			int len=strlen(tmp);
			if(strncmp(&tmp[len-3], " x", 2)==0) posprop=props;
			if(strncmp(&tmp[len-3], "nx", 2)==0) normprop=props;
			props++;
		}
		if(strcmp(tmp, "end_header\n")==0) break;
	}
	if(posprop==-1) {
		printf("ERROR: MeshObject::LoadPLY()- No vertex positions found\n");
		fclose(f);
		return false;
	}
	// Read verts
	int i=0;
	if(numverts>0) {
		numVertices = numverts;
		vec3* Vertexes = new vec3[numVertices];
		vec3* Normals = new vec3[numVertices];
		for(i=0; i< numVertices; i++) {
			fgets(tmp, 256, f);
			char *pch=strtok(tmp, " ");
			int prop=0;
			while(pch) {
				if(prop==posprop){ Vertexes[i].x=float(atof(pch)); 
					if(Vertexes[i].x < minX) minX = Vertexes[i].x; 
					if(Vertexes[i].x > maxX) maxX = Vertexes[i].x;
				}
				if(prop==posprop+1){ Vertexes[i].y=float(atof(pch));
					if(Vertexes[i].y < minY) minY = Vertexes[i].y;
					if(Vertexes[i].y > maxY) maxY = Vertexes[i].y;
				}
				if(prop==posprop+2){ Vertexes[i].z=float(atof(pch));
					if(Vertexes[i].z < minZ) minZ = Vertexes[i].z;
					if(Vertexes[i].z > maxZ) maxZ = Vertexes[i].z;
				}
				if(prop==normprop) Normals[i].x=float(atof(pch));
				if(prop==normprop+1) Normals[i].y=float(atof(pch));
				if(prop==normprop+2) Normals[i].z=float(atof(pch));
				pch=strtok(0, " ");
				prop++;
			}
		}
		// Now assign the arrays
		positions.assign(Vertexes, Vertexes + numVertices);
		normals.assign(Normals, Normals + numVertices);
		delete[] Vertexes;
		delete[] Normals;
	}
	// Read tris
	if(numtris>0) {
		//if(mtl==0) mtl=new LambertMaterial;
		int * Triangles = new int[3 * numtris];
		for(i=0; i<numtris; i++) {
			int count, i0, i1, i2;
			fscanf(f, "%d %d %d %d\n", &count, &i0, &i1, &i2);
			if(count!=3) {
				printf("ERROR: MeshObject::LoadPLY()- Only triangles are supported\n");
					fclose(f);
				return false;
			}
			Triangles[3 * i] = i0;
			Triangles[3 * i + 1] = i1;
			Triangles[3 * i + 2] = i2;
			triangles.push_back(i);
		}
		indices.assign(Triangles, Triangles + 3 * numtris);
		delete[] Triangles;
	}
	// Smooth
	if(normprop<0) smooth();
	// Close file
	fclose(f);
	cout << "Vertices count: " << numVertices << endl;
	cout << "Loaded " << numtris << " triangles from file " << filename << "!\n";
	return true;
}

void PLY::smooth(){
	cout << "Smoothing object...";
	for(int i=0; i< numVertices; i++)
		normals[i] = vec3(0);

	for(unsigned i=0; i< indices.size(); i += 3) {
		vec3 e1 = positions[indices[i + 1]] - positions[indices[i]];
		vec3 e2 = positions[indices[i + 2]] - positions[indices[i]];
		vec3 cr = cross(e1, e2);
		for(int j=0; j<3; j++)
			normals[indices[i + j]] += cr;
	}

	for(int i=0; i<numVertices; i++)
		normals[i] = normalize(normals[i]);
	cout << "done!" << endl;
}

void BoxTreeObject::Construct(PLY* obj){
	root = new TreeNode(obj, this);
	root->min = vec3(obj->minX, obj->minY, obj->minZ);
	root->max = vec3(obj->maxX, obj->maxY, obj->maxZ);
	root->Construct(obj, this, obj->triangles);
}

void TreeNode::Construct(PLY* obj, BoxTreeObject* tr, vector<int> indices){
	vector<int> leftV, rightV;
	int dim;
	if((max[0] - min[0]) > (max[1] - min[1])){ // X > Y
		if((max[0] - min[0]) > (max[2] - min[2])){ // X > Z
			dim = 0;
		}
		else{ // Z > X
			dim = 2;
		}
	}
	else{ // Y > X
		if((max[1] - min[1]) > (max[2] - min[2])){ // Y > Z
			dim = 1;
		}
		else{ // Z > Y
			dim = 2;
		}
	}
	left = new TreeNode(obj, tr);
	right = new TreeNode(obj, tr);
	float minLX = DEFAULT_LARGE_NUM, minLY = DEFAULT_LARGE_NUM, minLZ = DEFAULT_LARGE_NUM,
		  maxLX = -DEFAULT_LARGE_NUM, maxLY = -DEFAULT_LARGE_NUM, maxLZ = -DEFAULT_LARGE_NUM;
	float minRX = DEFAULT_LARGE_NUM, minRY = DEFAULT_LARGE_NUM, minRZ = DEFAULT_LARGE_NUM,
		  maxRX = -DEFAULT_LARGE_NUM, maxRY = -DEFAULT_LARGE_NUM, maxRZ = -DEFAULT_LARGE_NUM;
	float splitPoint = (max[dim] + min[dim]) / 2.f;
	vector<vec3>& pos = obj->getPositions();
	vector<int>& ind = obj->getIndices();
	for(unsigned int i = 0; i < indices.size(); i++){
		int offset = indices[i] * 3;
		vec3 center = (pos[ind[offset]] + pos[ind[offset + 1]] 
			+ pos[ind[offset + 2]]) / 3.f;
		if(center[dim] > splitPoint){ // Greater than split, put in right subtree
			rightV.push_back(indices[i]);
			if(pos[ind[offset]].x < minRX) minRX = pos[ind[offset]].x;
			if(pos[ind[offset]].y < minRY) minRY = pos[ind[offset]].y;
			if(pos[ind[offset]].z < minRZ) minRZ = pos[ind[offset]].z;
			if(pos[ind[offset]].x > maxRX) maxRX = pos[ind[offset]].x;
			if(pos[ind[offset]].y > maxRY) maxRY = pos[ind[offset]].y;
			if(pos[ind[offset]].z > maxRZ) maxRZ = pos[ind[offset]].z;

			if(pos[ind[offset + 1]].x < minRX) minRX = pos[ind[offset + 1]].x;
			if(pos[ind[offset + 1]].y < minRY) minRY = pos[ind[offset + 1]].y;
			if(pos[ind[offset + 1]].z < minRZ) minRZ = pos[ind[offset + 1]].z;
			if(pos[ind[offset + 1]].x > maxRX) maxRX = pos[ind[offset + 1]].x;
			if(pos[ind[offset + 1]].y > maxRY) maxRY = pos[ind[offset + 1]].y;
			if(pos[ind[offset + 1]].z > maxRZ) maxRZ = pos[ind[offset + 1]].z;

			if(pos[ind[offset + 2]].x < minRX) minRX = pos[ind[offset + 2]].x;
			if(pos[ind[offset + 2]].y < minRY) minRY = pos[ind[offset + 2]].y;
			if(pos[ind[offset + 2]].z < minRZ) minRZ = pos[ind[offset + 2]].z;
			if(pos[ind[offset + 2]].x > maxRX) maxRX = pos[ind[offset + 2]].x;
			if(pos[ind[offset + 2]].y > maxRY) maxRY = pos[ind[offset + 2]].y;
			if(pos[ind[offset + 2]].z > maxRZ) maxRZ = pos[ind[offset + 2]].z;
		}
		else{ // Smaller than split, put in left subtree
			leftV.push_back(indices[i]);
			if(pos[ind[offset]].x < minLX) minLX = pos[ind[offset]].x;
			if(pos[ind[offset]].y < minLY) minLY = pos[ind[offset]].y;
			if(pos[ind[offset]].z < minLZ) minLZ = pos[ind[offset]].z;
			if(pos[ind[offset]].x > maxLX) maxLX = pos[ind[offset]].x;
			if(pos[ind[offset]].y > maxLY) maxLY = pos[ind[offset]].y;
			if(pos[ind[offset]].z > maxLZ) maxLZ = pos[ind[offset]].z;

			if(pos[ind[offset + 1]].x < minLX) minLX = pos[ind[offset + 1]].x;
			if(pos[ind[offset + 1]].y < minLY) minLY = pos[ind[offset + 1]].y;
			if(pos[ind[offset + 1]].z < minLZ) minLZ = pos[ind[offset + 1]].z;
			if(pos[ind[offset + 1]].x > maxLX) maxLX = pos[ind[offset + 1]].x;
			if(pos[ind[offset + 1]].y > maxLY) maxLY = pos[ind[offset + 1]].y;
			if(pos[ind[offset + 1]].z > maxLZ) maxLZ = pos[ind[offset + 1]].z;

			if(pos[ind[offset + 2]].x < minLX) minLX = pos[ind[offset + 2]].x;
			if(pos[ind[offset + 2]].y < minLY) minLY = pos[ind[offset + 2]].y;
			if(pos[ind[offset + 2]].z < minLZ) minLZ = pos[ind[offset + 2]].z;
			if(pos[ind[offset + 2]].x > maxLX) maxLX = pos[ind[offset + 2]].x;
			if(pos[ind[offset + 2]].y > maxLY) maxLY = pos[ind[offset + 2]].y;
			if(pos[ind[offset + 2]].z > maxLZ) maxLZ = pos[ind[offset + 2]].z;
		}
	}
	left->min = vec3(minLX, minLY, minLZ);
	left->max = vec3(maxLX, maxLY, maxLZ);
	if(leftV.size() <= MAX_PER_BOX && leftV.size() > 0){
		left->triangleIndices.assign(leftV.data(), leftV.data() + leftV.size());
	}
	else if(leftV.size() > 0){
		left->Construct(obj, tr, leftV);
	}
	else if(leftV.size() == 0){
		int rInd = rightV[rightV.size() - 1];
		rightV.pop_back();
		left->triangleIndices.push_back(rInd);
		left->min.x = (pos[ind[rInd * 3]].x < pos[ind[rInd * 3 + 1]].x)?
			(pos[ind[rInd * 3]].x < pos[ind[rInd * 3 + 2]].x)?
			pos[ind[rInd * 3]].x : pos[ind[rInd * 3 + 2]].x :
			(pos[ind[rInd * 3 + 2]].x < pos[ind[rInd * 3 + 1]].x)?
			pos[ind[rInd * 3 + 2]].x : pos[ind[rInd * 3 + 1]].x;
		left->min.y = (pos[ind[rInd * 3]].y < pos[ind[rInd * 3 + 1]].y)?
			(pos[ind[rInd * 3]].y < pos[ind[rInd * 3 + 2]].y)?
			pos[ind[rInd * 3]].y : pos[ind[rInd * 3 + 2]].y :
			(pos[ind[rInd * 3 + 2]].y < pos[ind[rInd * 3 + 1]].y)?
			pos[ind[rInd * 3 + 2]].y : pos[ind[rInd * 3 + 1]].y;
		left->min.z = (pos[ind[rInd * 3]].z < pos[ind[rInd * 3 + 1]].z)?
			(pos[ind[rInd * 3]].z < pos[ind[rInd * 3 + 2]].z)?
			pos[ind[rInd * 3]].z : pos[ind[rInd * 3 + 2]].z :
			(pos[ind[rInd * 3 + 2]].z < pos[ind[rInd * 3 + 1]].z)?
			pos[ind[rInd * 3 + 2]].z : pos[ind[rInd * 3 + 1]].z;
		left->max.x = (pos[ind[rInd * 3]].x > pos[ind[rInd * 3 + 1]].x)?
			(pos[ind[rInd * 3]].x > pos[ind[rInd * 3 + 2]].x)?
			pos[ind[rInd * 3]].x : pos[ind[rInd * 3 + 2]].x :
			(pos[ind[rInd * 3 + 2]].x > pos[ind[rInd * 3 + 1]].x)?
			pos[ind[rInd * 3 + 2]].x : pos[ind[rInd * 3 + 1]].x;
		left->max.y = (pos[ind[rInd * 3]].y > pos[ind[rInd * 3 + 1]].y)?
			(pos[ind[rInd * 3]].y > pos[ind[rInd * 3 + 2]].y)?
			pos[ind[rInd * 3]].y : pos[ind[rInd * 3 + 2]].y :
			(pos[ind[rInd * 3 + 2]].y > pos[ind[rInd * 3 + 1]].y)?
			pos[ind[rInd * 3 + 2]].y : pos[ind[rInd * 3 + 1]].y;
		left->max.z = (pos[ind[rInd * 3]].z > pos[ind[rInd * 3 + 1]].z)?
			(pos[ind[rInd * 3]].z > pos[ind[rInd * 3 + 2]].z)?
			pos[ind[rInd * 3]].z : pos[ind[rInd * 3 + 2]].z :
			(pos[ind[rInd * 3 + 2]].z > pos[ind[rInd * 3 + 1]].z)?
			pos[ind[rInd * 3 + 2]].z : pos[ind[rInd * 3 + 1]].z;
	}
	right->min = vec3(minRX, minRY, minRZ);
	right->max = vec3(maxRX, maxRY, maxRZ);
	if(rightV.size() <= MAX_PER_BOX && rightV.size() > 0){
		right->triangleIndices.assign(rightV.data(), rightV.data() + rightV.size());
	}
	else if(rightV.size() > 0){
		right->Construct(obj, tr, rightV);
	}
	else if(rightV.size() == 0){
		int rInd = leftV[leftV.size() - 1];
		leftV.pop_back();
		right->triangleIndices.push_back(rInd);
		right->min.x = (pos[ind[rInd * 3]].x < pos[ind[rInd * 3 + 1]].x)?
			(pos[ind[rInd * 3]].x < pos[ind[rInd * 3 + 2]].x)?
			pos[ind[rInd * 3]].x : pos[ind[rInd * 3 + 2]].x :
			(pos[ind[rInd * 3 + 2]].x < pos[ind[rInd * 3 + 1]].x)?
			pos[ind[rInd * 3 + 2]].x : pos[ind[rInd * 3 + 1]].x;
		right->min.y = (pos[ind[rInd * 3]].y < pos[ind[rInd * 3 + 1]].y)?
			(pos[ind[rInd * 3]].y < pos[ind[rInd * 3 + 2]].y)?
			pos[ind[rInd * 3]].y : pos[ind[rInd * 3 + 2]].y :
			(pos[ind[rInd * 3 + 2]].y < pos[ind[rInd * 3 + 1]].y)?
			pos[ind[rInd * 3 + 2]].y : pos[ind[rInd * 3 + 1]].y;
		right->min.z = (pos[ind[rInd * 3]].z < pos[ind[rInd * 3 + 1]].z)?
			(pos[ind[rInd * 3]].z < pos[ind[rInd * 3 + 2]].z)?
			pos[ind[rInd * 3]].z : pos[ind[rInd * 3 + 2]].z :
			(pos[ind[rInd * 3 + 2]].z < pos[ind[rInd * 3 + 1]].z)?
			pos[ind[rInd * 3 + 2]].z : pos[ind[rInd * 3 + 1]].z;
		right->max.x = (pos[ind[rInd * 3]].x > pos[ind[rInd * 3 + 1]].x)?
			(pos[ind[rInd * 3]].x > pos[ind[rInd * 3 + 2]].x)?
			pos[ind[rInd * 3]].x : pos[ind[rInd * 3 + 2]].x :
			(pos[ind[rInd * 3 + 2]].x > pos[ind[rInd * 3 + 1]].x)?
			pos[ind[rInd * 3 + 2]].x : pos[ind[rInd * 3 + 1]].x;
		right->max.y = (pos[ind[rInd * 3]].y > pos[ind[rInd * 3 + 1]].y)?
			(pos[ind[rInd * 3]].y > pos[ind[rInd * 3 + 2]].y)?
			pos[ind[rInd * 3]].y : pos[ind[rInd * 3 + 2]].y :
			(pos[ind[rInd * 3 + 2]].y > pos[ind[rInd * 3 + 1]].y)?
			pos[ind[rInd * 3 + 2]].y : pos[ind[rInd * 3 + 1]].y;
		right->max.z = (pos[ind[rInd * 3]].z > pos[ind[rInd * 3 + 1]].z)?
			(pos[ind[rInd * 3]].z > pos[ind[rInd * 3 + 2]].z)?
			pos[ind[rInd * 3]].z : pos[ind[rInd * 3 + 2]].z :
			(pos[ind[rInd * 3 + 2]].z > pos[ind[rInd * 3 + 1]].z)?
			pos[ind[rInd * 3 + 2]].z : pos[ind[rInd * 3 + 1]].z;
	}
}

vec4 BoxTreeObject::intersection(vec3 c, Ray r, float& t, int d){
	if(root->intersect(r) < 0)
		return vec4(c, 1);
	else{
		float lt = -1, rt = -1;
		if(root->left)
			lt = root->left->intersect(r);
		if(root->right)
			rt = root->right->intersect(r);
		// Both negative
		if(lt < 0 && rt < 0){
			return vec4(c, 1);
		}
		// Only right positive
		else if(lt < 0){
			if(rt > 0 && rt < t){
				return root->right->intersection(c, r, t, d);
			}
			else
				return vec4(c, 1);
		}
		// Only left positive
		else if(rt < 0){
			if(lt > 0 && lt < t){
				return root->left->intersection(c, r, t, d);
			}
			else
				return vec4(c, 1);
		}
		// Both positive
		else{
			vec3 result = c;
			if(lt < rt){
				if(lt < t)
					result = vec3(root->left->intersection(result, r, t, d));
				else
					return vec4(c, 1);
				result = vec3(root->right->intersection(result, r, t, d));
			}
			else{
				if(rt < t)
					result = vec3(root->right->intersection(result, r, t, d));
				else
					return vec4(c, 1);
				result = vec3(root->left->intersection(result, r, t, d));
			}
			return vec4(result, 1);
		}
	}
}

float TreeNode::intersect(Ray r){
	// Purely empty node, do nothing
	if(!(triangleIndices.size() || left || right)) return -1;

	vec3 tMin = vec3(tree->compositeModel() * vec4(min, 1));
	vec3 tMax = vec3(tree->compositeModel() * vec4(max, 1));

	float tx1 = (tMin[0] - r.origin[0]) / r.direction[0];
	float tx2 = (tMax[0] - r.origin[0]) / r.direction[0];
	float ty1 = (tMin[1] - r.origin[1]) / r.direction[1];
	float ty2 = (tMax[1] - r.origin[1]) / r.direction[1];
	float tz1 = (tMin[2] - r.origin[2]) / r.direction[2];
	float tz2 = (tMax[2] - r.origin[2]) / r.direction[2];

	float tmin = glm::max(glm::max(glm::min(tx1, tx2), glm::min(ty1, ty2)), glm::min(tz1, tz2));
	float tmax = glm::min(glm::min(glm::max(tx1, tx2), glm::max(ty1, ty2)), glm::max(tz1, tz2));

	if(tmin <= tmax){
		if(tmin > 0){
			return tmin;
		}
		else if(tmax > 0){
			return tmax;
		}
	}
	return -1;
}

vec4 TreeNode::intersection(vec3& bg, Ray r, float& t, int d){
	// Check if leaf node. If it is, check per triangle
	if(triangleIndices.size()){
		vec3 returnCol = bg;
		mat4 compModel = tree->compositeModel();
		for(unsigned int i = 0; i < triangleIndices.size(); i++){
			vec3 a = vec3(compModel * vec4(o->vertexAt(o->indexAt(triangleIndices[i] * 3)), 1.f));
			vec3 b = vec3(compModel * vec4(o->vertexAt(o->indexAt(triangleIndices[i] * 3 + 1)), 1.f));
			vec3 c = vec3(compModel * vec4(o->vertexAt(o->indexAt(triangleIndices[i] * 3 + 2)), 1.f));
			vec3 normal = o->normalAt(o->indexAt(triangleIndices[i] * 3)) +
				o->normalAt(o->indexAt(triangleIndices[i] * 3 + 1)) +
				o->normalAt(o->indexAt(triangleIndices[i] * 3 + 2));
			normal = normalize(inverse(transpose(mat3(compModel))) * normal);
			returnCol = intersectTriangle(tree->getMaterial(), returnCol, r, a, b, c, t, normal, d);
		}
		return vec4(returnCol, 1.f);
	}
	// Not leaf, check children
	else{
		// For keeping track of time hit for left and right children
		float lt = -9, rt = -9;
		if(left)
			lt = left->intersect(r);
		if(right)
			rt = right->intersect(r);
		// Negative time for left, can't hit it
		if(lt < 0){
			// No hit if right also negative
			if(rt < 0)
				return vec4(bg, 1);
			// Hit occured, check the right subtree
			else{
				return right->intersection(bg, r, t, d);
			}
		}
		else{
			// Left positive, but right negative, so check left only
			if(rt < 0){
				return left->intersection(bg, r, t, d);
			}
			else{
				// Both hit, determine which to test first
				vec3 finalCol = bg;
				if(lt < rt){
					finalCol = left->intersection(bg, r, t, d);
					finalCol = right->intersection(finalCol, r, t, d);
				}
				else{
					finalCol = right->intersection(bg, r, t, d);
					finalCol = left->intersection(finalCol, r, t, d);
				}
				return vec4(finalCol, 1);
			}
		}
	}
}

bool TreeNode::intersects(Ray& r, vec4& max){
	// Check if leaf node. If it is, check per triangle
	if(triangleIndices.size()){
		mat4 compModel = tree->compositeModel();
		bool result = false;
		for(unsigned int i = 0; i < triangleIndices.size(); i++){
			vec3 a = vec3(compModel * vec4(o->vertexAt(o->indexAt(triangleIndices[i] * 3)), 1.f));
			vec3 b = vec3(compModel * vec4(o->vertexAt(o->indexAt(triangleIndices[i] * 3 + 1)), 1.f));
			vec3 c = vec3(compModel * vec4(o->vertexAt(o->indexAt(triangleIndices[i] * 3 + 2)), 1.f));
			result = intersectsTriangle(r, a, b, c, max);
			if(result)
				break;
		}
		return result;
	}
	// Not leaf, check children
	else{
		// For keeping track of time hit for left and right children
		float lt = -9, rt = -9;
		if(left)
			lt = left->intersect(r);
		if(right)
			rt = right->intersect(r);
		// Negative time for left, can't hit it
		if(lt < 0){
			// No hit if right also negative
			if(rt < 0)
				return false;
			// Hit occured, check the right subtree
			else{
				return right->intersects(r, max);
			}
		}
		else{
			// Left positive, but right negative, so check left only
			if(rt < 0){
				return left->intersects(r, max);
			}
			else{
				// Both hit, determine which to test first
				if(lt < rt){
					return left->intersects(r, max) || right->intersects(r, max);
				}
				else{
					return right->intersects(r, max) || left->intersects(r, max);
				}
			}
		}
	}
}

bool BoxTreeObject::intersects(Ray& r, vec4& max){
	if(root->intersect(r) < 0)
		return false;
	else{
		float lt = -1, rt = -1;
		if(root->left)
			lt = root->left->intersect(r);
		if(root->right)
			rt = root->right->intersect(r);
		// Both negative
		if(lt < 0 && rt < 0){
			return false;
		}
		// Only right positive
		else if(lt < 0){
			if(rt > 0){
				return root->right->intersects(r, max);
			}
			else
				return false;
		}
		// Only left positive
		else if(rt < 0){
			if(lt > 0){
				return root->left->intersects(r, max);
			}
			else
				return false;
		}
		// Both positive
		else{
			if(lt < rt){
				return root->left->intersects(r, max) || root->right->intersects(r, max);
			}
			else{
				return root->right->intersects(r, max) || root->left->intersects(r, max);
			}
		}
	}
}