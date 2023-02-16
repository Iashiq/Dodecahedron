//=============================================================================================
// Computer Graphics Second Homework
// I hereby declare that the homework has been made by me, including the problem interpretation,
// algorithm selection, and coding. Should I use materials and programs not from the course webpage, 
// the sources are clearly indentified as a comment in the code. 
//=============================================================================================
#include "framework.h"

// vertex shader
const char* vertexSource = R"(
	#version 330
    precision highp float;
	
	uniform vec3 wLookAt, wRight, wUp; 
	layout(location = 0) in vec2 cCamWindowVertex;	
	out vec3 p;

	void main() 
	{
		gl_Position = vec4(cCamWindowVertex, 0, 1);
		p = wLookAt + wRight * cCamWindowVertex.x + wUp * cCamWindowVertex.y;
	}
)";
//fragment shader
const char* fragmentSource = R"(
	#version 330
    precision highp float;

	uniform vec3 kdo[2], kel[2], F0;
	uniform vec3 cameraEye, vertices[20];
    const vec3 ka = vec3(0.4f, 0.2f, 0.5f);
    const vec3 positionOfLight = vec3(0.2f, 0.1f, 0.12f);
	const vec3 background = vec3(0.1f, 0.4f, 0.5f);/////////////////03.,0,4,0,5/
	const vec3 color = vec3(0.3f, 0.9f, 0.9f);
	
	
	
	
    const int noOfEdges = 60;
	

	const int totalFaces = 12;
    uniform int noOfPlanes[totalFaces * 36];
	uniform int top;
    const float shininess = 150.0f;
    const float epsilon = 0.01f;///0.7f
    const int depth = 5;



	struct Hit
	{
		float t;
		vec3 position, normalVector;
		int material;
	};

	struct Ray
	{
		vec3 originOfLight, directionOfLight, weight;
	};



	Hit drawEllipsoid(Ray ray, Hit hit, vec3 center, float a0, float b0, float c0, float r)
	{
		vec3 distance = ray.originOfLight - center;

    ////////////////To Calculate the value of a,b and c, I got help from the source below
    ///////////////https://github.com/roger-vertices/raytracing/blob/master/ellipsoid.cpp

		float a =  (ray.directionOfLight.x * ray.directionOfLight.x / a0 + ray.directionOfLight.y * ray.directionOfLight.y / b0 + ray.directionOfLight.z * ray.directionOfLight.z / c0);
		float b = 2.0f * (distance.x * ray.directionOfLight.x / a0 +  distance.y * ray.directionOfLight.y / b0 +  distance.z * ray.directionOfLight.z / c0); 
		float c =  (distance.x * distance.x / a0 + distance.y * distance.y / b0 + distance.z * distance.z / c0) - (r * r);               

		float discriminent = sqrt((b * b) - (4.0f * a * c));
		float root1 = (-b + discriminent) / 2.0f / a;
		float root2 = (-b - discriminent) / 2.0f / a;
		if (root1 <= 0)
		{
			return hit;
		}

		hit.t = (root2 > 0) ? root2 : root1;
		vec3 position = ray.originOfLight + ray.directionOfLight * hit.t;
		if (position.y <= 0.5)
		{
			hit.position = position;
			hit.normalVector = normalize(vec3((hit.position.x - center.x) * 2 / a0, (hit.position.y - center.y) * 2 / b0,(hit.position.z - center.z) * 2 / c0));
			hit.material = 99;
			return hit;
		}
		else 
		{
			hit.t = -1;
			return hit;
		}
	}

	void connectPlane(int i, float width, out vec3 p, out vec3 normalVector)
	{
		vec3 point1 = vertices[noOfPlanes[3*i]-1];
		vec3 point2 = vertices[noOfPlanes[3*i+1]-1];
        vec3 point3 = vertices[noOfPlanes[3*i+2]-1];

		normalVector = cross(point2-point1, point3-point1);
		if(dot(point1, normalVector) < 0)
			normalVector = -normalVector;
		
		p = point1 * width + vec3(0,0,0.001f);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////To draw dodecahedron, I got little help from the lecture video/////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Hit drawDodecahedron(Ray ray, Hit hit, float width, int material)
	{
		vec3 a, b;
		int arrayOfEdges[noOfEdges] = int[] (1, 2, 16, 5, 13, 1, 13, 9, 10, 14, 1, 14, 6, 15, 2, 2, 15, 11, 12, 16, 3, 4, 18, 8, 17, 3, 17, 12, 11, 20,3, 20, 7, 19, 4, 19, 10, 9, 
            18,  4, 16, 12, 17, 8, 5,  5, 8, 18, 9, 13, 14, 10, 19, 7, 6, 6, 7, 20, 11, 15);

		for(int i = 0; i < totalFaces; i++) 
		{
			vec3 point1, normalVector;
			connectPlane(i, width, point1, normalVector);
			float root1=abs(dot(normalVector, ray.directionOfLight)) > epsilon ? dot(point1 - ray.originOfLight, normalVector) / dot(normalVector, ray.directionOfLight) : -1;
			if (root1 <= epsilon || (root1 > hit.t && hit.t > 0)) continue;
			vec3 pointOfIntersection = ray.originOfLight + ray.directionOfLight * root1;
			bool isOut = false;

			for(int j = 0; j< totalFaces; j++) 
			{
				if (i == j) continue;
				vec3 point2, n;
				connectPlane(j, width, point2, n);
				if (dot(n, pointOfIntersection - point2) > 0)
				{
					isOut = true;
					break;
				}

			}

			if (!isOut) 
			{
				hit.t = root1;
				hit.position = pointOfIntersection;
				hit.normalVector = normalize(normalVector);
				hit.material = 99;  ///99//0////////////////////////////////////////////////////////
				int l = 0;
                while(l<depth)
				{
					a = vertices[arrayOfEdges[i * 5 + l] - 1];
					b = vertices[arrayOfEdges[i * 5 + (l == 4 ? 0 : l + 1)] - 1];
					
					vec3 P = ray.originOfLight + ray.directionOfLight * hit.t;
					float k = dot(P - a, normalize(b - a));
					vec3 t = P - a - normalize(b - a) * k;
					float dist = dot(P - a, normalize(t));
					if(dist < 0.18f)
					{
						hit.material = material;
						break;
					}
                   l++;
				}
			}
		}
		return hit;
	}


	Hit drawOrbifold(Ray ray) 
	{
		Hit hit;
		hit.t = -0.5;
		hit = drawEllipsoid(ray, hit, vec3(0.01f, 0.07f, 0.04f), 1.5f, 2.1f, 0.99f, 0.09f);
		hit = drawDodecahedron(ray, hit, 0.95f, 1);
		if (dot (ray.directionOfLight, hit.normalVector) > 0) hit.normalVector = hit.normalVector * (-1);
		return hit;
	}


	vec3 rayTrace(Ray ray)
	{
		vec3 radiance = vec3(0, 0, 0);
        int i=0;
		while(i<5)
		{
			Hit hit = drawOrbifold(ray);
			if (hit.t < 0) break;
			if (hit.material < 2) 
			{
				vec3 lightdir = normalize(positionOfLight - hit.position);
				float cosTheta = dot(hit.normalVector, lightdir);
				if (cosTheta > 0) 
				{
					vec3 LeIn = color / dot(positionOfLight - hit.position, positionOfLight - hit.position);
					radiance += ray.weight * LeIn * kdo[hit.material] * cosTheta;
					vec3 halfway = normalize(-ray.directionOfLight + lightdir);
					float cosDelta = dot(hit.normalVector, halfway);
					if (cosDelta > 0) radiance += ray.weight * LeIn * kel[hit.material] * pow (cosDelta, shininess);
				}
					ray.weight *= ka;
					break;
			}

			if(hit.material == 99)
			{
				ray.weight *= F0 + (vec3(1, 1, 1) - F0) * pow( dot(-ray.directionOfLight, hit.normalVector), 5);
				ray.originOfLight = hit.position + hit.normalVector * epsilon;
				ray.directionOfLight = reflect(ray.directionOfLight, hit.normalVector);
				i++;
			}
			else if(hit.material == 0)
			{
				ray.originOfLight = hit.position + hit.normalVector * epsilon;
				ray.directionOfLight = reflect(ray.directionOfLight, hit.normalVector);
				
			}
			i++;
		}
		radiance += ray.weight*background;
		return radiance;
	}

	in vec3 p;
	out vec4 fragmentColor;

	void main()
	{
		Ray ray;
		ray.originOfLight = cameraEye;
		ray.directionOfLight = normalize(p-cameraEye);
		ray.weight = vec3(4,2,2); ///(1 1 1)//2,3,3,// 4,2,2
		fragmentColor = vec4(rayTrace(ray), 1);
	}



)";


struct Camera
{
	vec3 eye, lookat, right, pvup, rvup;
	float fov = 110 * (float)M_PI / 180;

	Camera() : eye(0.3, 0.2, 0.5), pvup(0, 0, 2), lookat(0, 0, 0) { set(); }

	void set()
	{
		vec3 w = eye - lookat;
		float f = length(w);
		right = normalize(cross(pvup, w)) * f * tanf(fov / 2);
		rvup = normalize(cross(w, right)) * f * tanf(fov / 2);
	}

	void animate(float t)
	{
		float r = sqrtf(eye.x * eye.x + eye.y * eye.y);
		eye = vec3(r * cos(t) + lookat.x, r * sin(t) + lookat.y, lookat.z);
		set();
	}

	void Step(float step)
	{
		eye = normalize(eye + pvup * step) * length(eye);
		set();
	}
};

GPUProgram shader;
Camera camera;
bool animate = true;



float F(float n, float k)
{
	return((n - 1) * (n - 1) + k * k) / ((n + 1) * (n + 1) + k * k);
}

// Initialization, create an OpenGL context
void onInitialization() {
	glViewport(0, 0, windowWidth, windowHeight);


	unsigned int vao, vbo;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);

	float vertexCoords[] = { -1, -1, 1, -1, 1, 1, -1, 1 };
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexCoords), vertexCoords, GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);

	shader.create(vertexSource, fragmentSource, "fragmentColor");
	shader.setUniform(1, "top");

	const float g = 0.618f;
	const float G = 1.618f;

	std::vector<vec3> vertices = {   vec3(0, 0.618, 1.618),	vec3(0, -0.618, 1.618), vec3(0, -0.618, -1.618), vec3(0, 0.618, -1.618), vec3(1.618, 0, 0.618),	vec3(-1.618, 0, 0.618),	vec3(-1.618, 0, -0.618), vec3(1.618, 0, -0.618),
		vec3(0.618, 1.618, 0), vec3(-0.618, 1.618, 0), vec3(-0.618, -1.618, 0), vec3(0.618, -1.618, 0), vec3(1, 1, 1), vec3(-1, 1, 1), vec3(-1, -1, 1), vec3(1, -1, 1), vec3(1, -1, -1), vec3(1, 1, -1),
		vec3(-1, 1, -1),vec3(-1, -1, -1) };


	for (int i = 0; i < vertices.size(); i++)
	{
		shader.setUniform(vertices[i], "vertices[" + std::to_string(i) + "]");
	}

	std::vector<int> noOfPlanes = {   1,   2,  16, 1,  13,   9, 1,  14,   6, 2,  15,  11, 3,   4,  18, 3,  17,  12, 
	                           	3,  20,   7, 19, 10,   9, 16, 12,  17, 5,   8,  18, 14, 10,  19, 6,   7,  20 };


	for (int i = 0; i < noOfPlanes.size(); i++)
	{
		shader.setUniform(noOfPlanes[i], "noOfPlanes[" + std::to_string(i) + "]");
	}
	shader.setUniform(vec3(0.1f, 0.2f, 0.4f), "kdo[0]");
	shader.setUniform(vec3(1.5, 0.6f, 0.4f), "kdo[1]");
	shader.setUniform(vec3(5, 5, 5), "kel[0]");
	shader.setUniform(vec3(1, 1, 1), "kel[1]");
	shader.setUniform(vec3(F(0.17, 3.1), F(0.35, 2.7), F(1.5, 1.9)), "F0");
}

// Window has become invalid: Redraw
void onDisplay()
{
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	shader.setUniform(camera.eye, "cameraEye");
	shader.setUniform(camera.lookat, "wLookAt");
	shader.setUniform(camera.right, "wRight");
	shader.setUniform(camera.rvup, "wUp");
	glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

	//...fullScreenTexturedQuad->Draw();
	glutSwapBuffers();									// exchange the two buffers
}

// Key of ASCII code pressed
void onKeyboard(unsigned char key, int pX, int pY) {
}

// Key of ASCII code released
void onKeyboardUp(unsigned char key, int pX, int pY) {

}

// Mouse click event
void onMouse(int button, int state, int pX, int pY) {
}

// Move mouse with key pressed
void onMouseMotion(int pX, int pY) {
}

// Idle event indicating that some time elapsed: do animation here
void onIdle()
{
	if (animate) camera.animate(glutGet(GLUT_ELAPSED_TIME) / 1000.0f);
	glutPostRedisplay();
}
