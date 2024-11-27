#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <TextureMap.h>
#include <RayTriangleIntersection.h>

#define __INT_MAX__ 2147483647
#define pi 3.14159265359
#define R 0.23
#define FOCAL_LEN 2.0
#define WIDTH 320
#define HEIGHT 240
#define OBJ_PATH "cornell-box.obj"

glm::vec3 camPos = glm::vec3(0.0f, 0.0f, 4.0f);
glm::mat3 camMat = glm::mat3(glm::vec3(1.0f, 0.0f, 0.0f), 
							glm::vec3(0.0f, 1.0f, 0.0f), 
							glm::vec3(0.0f, 0.0f, 1.0f));
glm::vec3 light = glm::vec3(0.0, 0.6, 0.0);
float deg = 0.01f;
bool orbit = false;
int mode = 3;

std::unordered_map<std::string, Colour> parseMtl(const std::string& filename) {
    std::unordered_map<std::string, Colour> mtls;
    std::ifstream inputStream(filename);
    std::string line;

    std::string currentColor;

    while (true) {
        getline(inputStream, line);
        if (inputStream.fail()) break;
        // split line
        std::vector<std::string> splittedLine = split(line, ' ');
        if (splittedLine[0] == "newmtl") {
            currentColor = splittedLine[1];
        }
        else if (splittedLine[0] == "Kd") {
            float red = std::stof(splittedLine[1]) * 255;
            float green = std::stof(splittedLine[2]) * 255;
            float blue = std::stof(splittedLine[3]) * 255;
            Colour c = Colour(currentColor, red, green, blue);
            mtls[currentColor] = c;
        }
    }

    return mtls;
}

std::vector<ModelTriangle> parseObj(const std::string& filename) {
	std::ifstream inputStream(filename);
	std::string line;

    std::vector<ModelTriangle> triangles, triangle;
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> faces;
    std::unordered_map<std::string, Colour> colours;
    std::string currentColor;
    std::string mtlFile;

    while (true) {
        getline(inputStream, line);
        if (inputStream.fail()) break;
        // split line
        std::vector<std::string> splittedLine = split(line, ' ');
        if (splittedLine[0] == "v") {
            float x = std::stof(splittedLine[1]) * 0.35;
            float y = std::stof(splittedLine[2]) * 0.35;
            float z = std::stof(splittedLine[3]) * 0.35;
            vertices.push_back(glm::vec3(x, y, z));
        }
        else if (splittedLine[0] == "f") {
            Colour triColour = colours[currentColor];
            glm::vec3 first = vertices[std::stoi(splittedLine[1]) - 1];
            glm::vec3 second = vertices[std::stoi(splittedLine[2]) - 1];
            glm::vec3 third = vertices[std::stoi(splittedLine[3]) - 1];
            glm::vec3 edge1 (first - second); 
            glm::vec3 edge2 (third - second);
            // normal vector to triangle
            glm::vec3 normal = glm::cross(edge1, edge2);
            ModelTriangle triangle (first, second, third, triColour);
            triangle.normal = normal;
            triangles.push_back(triangle);
            
        }
        else if (splittedLine[0] == "o") {
        }
        else if (splittedLine[0] == "usemtl") {
            currentColor = splittedLine[1];
        }
        else if (splittedLine[0] == "mtllib") {
            colours = parseMtl(splittedLine[1]);
        }
    }
    
    return triangles;
}

CanvasTriangle randTri(){
	CanvasPoint v1(rand() %319, rand() %239);
	CanvasPoint v2(rand() %319, rand() %239);
	CanvasPoint v3(rand() %319, rand() %239);
	CanvasTriangle tri(v1, v2, v3);

	return tri;
}

Colour randCol() {
    Colour col = {rand()%256, rand()%256, rand()%256};
    return col;
}

void resetWindow(DrawingWindow& window){
	window.clearPixels();
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	float space1 = (to.x - from.x) / (numberOfValues - 1);
	float space2 = (to.y - from.y) / (numberOfValues - 1);
	float space3 = (to.z - from.z) / (numberOfValues - 1);

	std::vector<glm::vec3> result;
   	for (int i = 0; i < numberOfValues; ++i)
    	result.push_back(glm::vec3 (from.x + i * space1, from.y + i * space2, from.z + i * space3 ));

	return result;
}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	float spacing = (to - from) / (numberOfValues - 1);
	std::vector<float> result;
   	for (int i = 0; i < numberOfValues; ++i)
    	result.push_back(from + i * spacing);

	return result;
}

std::vector<CanvasPoint> interpolateCanvasPointWithTexture(CanvasPoint from, CanvasPoint to, int numberOfValues) {
	std::vector<CanvasPoint> steps;
	float stepSizeX = (to.x - from.x) / (numberOfValues - 1);
	float stepSizeY = (to.y - from.y) / (numberOfValues - 1);
	float stepSizeTextureX = (to.texturePoint.x - from.texturePoint.x) / (numberOfValues - 1);
	float stepSizeTextureY = (to.texturePoint.y - from.texturePoint.y) / (numberOfValues - 1);

	for (int i = 0; i < numberOfValues; i++) {
		CanvasPoint n = CanvasPoint(from.x + i * stepSizeX, from.y + i * stepSizeY);
		n.texturePoint.x = from.texturePoint.x + i * stepSizeTextureX;
		n.texturePoint.y = from.texturePoint.y + i * stepSizeTextureY;
		steps.push_back(n);
	}
	return steps;
}

float calcBrightness(glm::vec3 light, glm::vec3 camPosition, glm::vec3 intersectionPoint, glm::vec3 normal) {
    float brightness = (1 / (4 * pi * R * R));
    glm::vec3 surfaceLight = light - intersectionPoint;
    float angle = glm::dot(normal, surfaceLight);

    //check value between 0-1
    float intensity = (brightness * angle);
    if (intensity > 1) {
        intensity = 1;
    } else if (intensity < 0.1) {
        intensity = 0.1;
    }

    return intensity;
}

RayTriangleIntersection getClosestValidIntersection(glm::vec3 camPosition, glm::vec3 rayDir, std::vector<ModelTriangle> modelTriangles) {
	RayTriangleIntersection minimumTintersection = RayTriangleIntersection(glm::vec3(),  2147483648, ModelTriangle(), 0);
	RayTriangleIntersection validIntersection;
	
	for (int i = 0; i < modelTriangles.size(); i++){
		ModelTriangle triangle = modelTriangles[i];

		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = camPosition - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDir, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

		//r = s + t * d
		glm::vec3 intersectionPoint = camPosition + rayDir * possibleSolution.x;
		// check if intersectionPoint inside triangle
		bool inU = 0 <= possibleSolution.y && possibleSolution.y <= 1;
		bool inV = 0 <= possibleSolution.z && possibleSolution.z <= 1;
		bool inUV = possibleSolution.y + possibleSolution.z <= 1;
		bool infrontCam = possibleSolution.x > 0;
		int validIntersectionPoint = (inU && inV && inUV && infrontCam) ? 1 : 0;
		RayTriangleIntersection validIntersection = RayTriangleIntersection(intersectionPoint, possibleSolution.x, triangle, validIntersectionPoint);

		if (validIntersection.triangleIndex == 1 && validIntersection.distanceFromCamera < minimumTintersection.distanceFromCamera){
			minimumTintersection = validIntersection;
			minimumTintersection.triangleIndex = i;
		}

	}

	return minimumTintersection;
}

bool shadowEffect(glm::vec3 currPoint, size_t currPointIndex, std::vector<ModelTriangle> modelTriangles) {
	glm::vec3 light = glm::vec3(0.0, 0.45, 0.5);
	glm::vec3 rayDir = glm::normalize(currPoint - light);
	RayTriangleIntersection shadowIntersection = getClosestValidIntersection(light, rayDir, modelTriangles);
	return shadowIntersection.triangleIndex == currPointIndex;
}

void drawRayTraceAngle(DrawingWindow& window, glm::vec3 camPosition, glm::mat3 camMat) {
	std::vector<ModelTriangle> modelTriangles = parseObj(OBJ_PATH);
	uint32_t colourNumeric = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(0);

	for (int j = 0; j < HEIGHT; j++) {
		for (int i = 0; i < WIDTH; i++) {

			// find the ray that coresponds with the point (i,j)
			glm::vec3 rayDir, rayDirN;
			glm::vec3 camSpace;
			int z = 1;
			int scale = 45;

			camSpace.x = z * (i - WIDTH / 2) / (scale * FOCAL_LEN);
			camSpace.y = -z * (j - HEIGHT / 2) / (scale * FOCAL_LEN);
			camSpace.z = 0;
			rayDirN = camSpace - camPosition;
			rayDir = glm::normalize(rayDirN * camMat);

			// find the intersection point
			RayTriangleIntersection intersectionPoint = getClosestValidIntersection(camPosition, rayDir, modelTriangles);

			// color the window at (j,i)
			Colour currCol = intersectionPoint.intersectedTriangle.colour;

			//look for shadow
			if (shadowEffect(intersectionPoint.intersectionPoint, intersectionPoint.triangleIndex, modelTriangles)) {
				colourNumeric = (255 << 24) + (int(currCol.red) << 16) + (int(currCol.green) << 8) + int(currCol.blue);
			}
			//proximity light calculation
			else {
				float brightness = calcBrightness(light, camPos, intersectionPoint.intersectionPoint, intersectionPoint.intersectedTriangle.normal);
				colourNumeric = (255 << 24) + (int(currCol.red * brightness) << 16) + (int(currCol.green * brightness) << 8) + int(currCol.blue * brightness);
			}

			window.setPixelColour(i, j, colourNumeric);
		}
	}
}

void drawRayTraceSoft(DrawingWindow& window, glm::vec3 camPosition, glm::mat3 camMat) {
	std::vector<ModelTriangle> modelTriangles = parseObj(OBJ_PATH);
	uint32_t colourNumeric = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(0);

	for (int j = 0; j < HEIGHT; j++) {
		for (int i = 0; i < WIDTH; i++) {

			// find the ray that coresponds with the point (i,j)
			glm::vec3 rayDir, rayDirN;
			glm::vec3 camSpace;
			int z = 1;
			int scale = 45;
			camSpace.x = z * (i - WIDTH / 2) / (scale * FOCAL_LEN);
			camSpace.y = -z * (j - HEIGHT / 2) / (scale * FOCAL_LEN);
			camSpace.z = 0;
			rayDirN = camSpace - camPosition;
			rayDir = glm::normalize(rayDirN * camMat);

			// find the intersection point
			RayTriangleIntersection intersectionPoint = getClosestValidIntersection(camPosition, rayDir, modelTriangles);

			// color the window at (j,i)
			Colour currCol = intersectionPoint.intersectedTriangle.colour;

			//look for shadow
			if (shadowEffect(intersectionPoint.intersectionPoint, intersectionPoint.triangleIndex, modelTriangles)) {
				colourNumeric = (255 << 24) + (int(currCol.red) << 16) + (int(currCol.green) << 8) + int(currCol.blue);
			}
			//proximity light calculation
			else {
				float brightness = glm::length(rayDirN) * (1 / (4 * pi * R * R));
				colourNumeric = (255 << 24) + (int(currCol.red * brightness) << 16) + (int(currCol.green * brightness) << 8) + int(currCol.blue * brightness);
			}

			window.setPixelColour(i, j, colourNumeric);
		}
	}
}

void drawRayTraceHard(DrawingWindow& window, glm::vec3 camPosition, glm::mat3 camMat) {
	std::vector<ModelTriangle> modelTriangles = parseObj(OBJ_PATH);
	uint32_t colourNumeric, black = (255 << 24) + (int(0) << 16) + (int(0) << 8) + int(0);

	for (int j = 0; j < HEIGHT; j++) {
		for (int i = 0; i < WIDTH; i++) {

			// find the ray that coresponds with the point (i,j)
			glm::vec3 rayDir;
			glm::vec3 camSpace;
			int z = 1;
			int scale = 45;
			camSpace.x = z * (i - WIDTH / 2) / (scale * FOCAL_LEN);
			camSpace.y = -z * (j - HEIGHT / 2) / (scale * FOCAL_LEN);
			camSpace.z = 0;
			rayDir = camSpace - camPosition;
			rayDir = glm::normalize(rayDir * camMat);

			// find the intersection point
			RayTriangleIntersection intersectionPoint = getClosestValidIntersection(camPosition, rayDir, modelTriangles);

			// color the window at (j,i)
			Colour currCol = intersectionPoint.intersectedTriangle.colour;

			//look for shadow
			if (shadowEffect(intersectionPoint.intersectionPoint, intersectionPoint.triangleIndex, modelTriangles)) {
				colourNumeric = (255 << 24) + (int(currCol.red) << 16) + (int(currCol.green) << 8) + int(currCol.blue);
			}
			//proximity light calculation
			else {
				colourNumeric = black;
			}

			window.setPixelColour(i, j, colourNumeric);
		}
	}
}

glm::vec3 moveUp(DrawingWindow& window, glm::vec3 camPos) {
	resetWindow(window);
	camPos.y = camPos.y - 1;
	return camPos;
}

glm::vec3 moveDown(DrawingWindow& window, glm::vec3 camPos) {
	resetWindow(window);
	camPos.y = camPos.y + 1;
	return camPos;
}

glm::mat3 rotateXaxis(DrawingWindow& window, glm::mat3 camMat) {
	resetWindow(window);
	float deg = 0.01;
	glm::mat3 xAxisMat = glm::mat3(glm::vec3(1.0f, 0.0f, 0.0f),
									glm::vec3(0.0f, cos(deg), sin(deg)),
									glm::vec3(0.0f, -sin(deg), cos(deg)));
	camMat = xAxisMat * camMat;
	return camMat;
}

glm::mat3 rotateYaxis(DrawingWindow& window, glm::mat3 camMat) {
	resetWindow(window);
	float deg = 0.01;
	glm::mat3 yAxisMat = glm::mat3(glm::vec3(cos(deg), 0.0f, sin(deg)),
									glm::vec3(0.0f, 1.0f, 0.0f),
									glm::vec3(-sin(deg), 0.0f, cos(deg)));
	camMat = yAxisMat * camMat;
	return camMat;
}

glm::mat3 lookAt(glm::vec3& camPos, glm::mat3& camMat) {
	glm::vec3 forward = glm::normalize(camPos - glm::vec3(0.0f,0.0f,0.0f));
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0.0f,1.0f,0.0f), forward));
	glm::vec3 up = glm::normalize(glm::cross(forward, right));

	camMat[0] = right;
	camMat[1] = up;
	camMat[2] = forward;
	return camMat;
}

void drawLine(CanvasPoint from, CanvasPoint to, DrawingWindow &window, Colour colour){
	float Xdiff = to.x - from.x;
	float Ydiff = to.y - from.y;
	float stepNum = fmax(abs(Xdiff),abs(Ydiff));
	float XStepSize = Xdiff/stepNum;
	float YStepSize = Ydiff/stepNum;
	for (float i = 0.0; i <= stepNum; i++){
		float x = from.x + (XStepSize * i);
		float y = from.y + (YStepSize * i);
		uint32_t colourRandom = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
		window.setPixelColour(round(x),round(y), colourRandom);
	}
}

void drawHorizontalTextureLine(CanvasPoint from, CanvasPoint to, TextureMap tm, DrawingWindow& window) {
	// std::cout << "!!New Line" << std::endl;
	std::vector<CanvasPoint> canvasLine = interpolateCanvasPointWithTexture(from, to, to.x - from.x);
	if (canvasLine.size() == 1) {
		uint32_t colourNumeric = tm.pixels[from.texturePoint.x + from.texturePoint.y * tm.width];
		window.setPixelColour(round(from.x), round(from.y), colourNumeric);
		return;
	}
	//std::cout << "Draw From: " << from << " t.x: " << from.texturePoint.x << " t.y: " << from.texturePoint.y << std::endl;
	//std::cout << "Draw To: " << to << " t.x: " << to.texturePoint.x << " t.y: " << to.texturePoint.y << std::endl;
	for (auto c: canvasLine) {
		// std::cout << "Draw: " << c << " t.x: " << c.texturePoint.x << " t.y: " << c.texturePoint.y << std::endl;
		int x = c.texturePoint.x;
		int y = c.texturePoint.y;
		uint32_t colourNumeric = tm.pixels.at(x + y * tm.width);
		//colourNumeric = (0xFF << 24) + (int(float(c.texturePoint.x) / 3.0) << 16) + (int(float(c.texturePoint.y) / 3.0) << 8);
		window.setPixelColour(round(c.x), round(c.y), colourNumeric);
	}
}

CanvasPoint fourthPoint(CanvasPoint top, CanvasPoint mid, CanvasPoint bottom, bool texture = false, bool depth = false) {
	if (mid.y == top.y) {
		return top;
	}
	else if (mid.y == bottom.y) {
		return bottom;
	}
	int fm = mid.y - top.y;
	int l0 = std::max(bottom.x, top.x) - std::min(bottom.x, top.x);
	int l2 = bottom.y - top.y;
	float xFourth, m0 = (fm * l0) / l2;
	if (top.x > bottom.x) {
		xFourth = top.x - m0;
	}
	else {
		xFourth = top.x + m0;
	}
	float yFourth = mid.y;
	CanvasPoint n = CanvasPoint(xFourth, yFourth);
	std::vector<CanvasPoint> fullSide;
	if (texture) {
		int numOfValues = bottom.y - top.y;
		fullSide = interpolateCanvasPointWithTexture(top, bottom, numOfValues);
		int newPointIndex = mid.y - top.y;
		n.texturePoint.x = fullSide[newPointIndex].texturePoint.x;
		n.texturePoint.y = fullSide[newPointIndex].texturePoint.y;
	}

	return n;
	}

CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 camPos, glm::vec3 vertexPosition, float focalLength, glm::mat3 camMat) {
	CanvasPoint canvasPoint = CanvasPoint();
	glm::vec3 cameraToVertex = vertexPosition - camPos;
	vertexPosition = cameraToVertex * camMat;
	canvasPoint.x = ((focalLength / vertexPosition.z) * -vertexPosition.x) * 160 + (WIDTH / 2);
	canvasPoint.y = ((focalLength / vertexPosition.z) * vertexPosition.y) * 160 + (HEIGHT / 2);
	canvasPoint.depth = -vertexPosition.z;
	return canvasPoint;
}	

void drawTextureTriangle(CanvasTriangle t, DrawingWindow& window, TextureMap tm) {
	CanvasPoint top = t[0];
	CanvasPoint mid = t[1];
	CanvasPoint bottom = t[2];

	// divide to 2 triangles (find 4th ver)
	CanvasPoint r, l;
	CanvasPoint n = fourthPoint(top, mid, bottom, true, false);
	
	// find right and left points
	if (n.x > mid.x) {
		r = n;
		l = mid;
	}
	else {
		r = mid;
		l = n;
	}
	
	// color top triangle
	float topTriangleSectionHeight = r.y - top.y;
	std::vector<CanvasPoint> rightSide = interpolateCanvasPointWithTexture(top, r, topTriangleSectionHeight);
	std::vector<CanvasPoint> leftSide = interpolateCanvasPointWithTexture(top, l, topTriangleSectionHeight);
	for (int i = 0; i < topTriangleSectionHeight; i++) {
		CanvasPoint f = leftSide.at(i);
		CanvasPoint t = rightSide.at(i);
		drawHorizontalTextureLine(f, t, tm, window);
	}
	// color bottom triangle
	rightSide = interpolateCanvasPointWithTexture(r, bottom, bottom.y - r.y);
	leftSide = interpolateCanvasPointWithTexture(l, bottom, bottom.y - l.y);
	for (int i = mid.y; i < bottom.y; i++) {
		CanvasPoint f = leftSide.at(i - mid.y);
		CanvasPoint t = rightSide.at(i - mid.y);
		drawHorizontalTextureLine(f, t, tm, window);
	}
}

void drawTriangle(CanvasTriangle triangle, DrawingWindow &window, Colour colour){
	drawLine(triangle[0], triangle[1], window, colour);
	drawLine(triangle[1], triangle[2], window, colour);
	drawLine(triangle[2], triangle[0], window, colour);

}

void drawFilTri(CanvasTriangle t, DrawingWindow& window, Colour colour) {
	// sort top to bottom
	CanvasPoint a[] = { t.v0(), t.v1(), t.v2() };
	std::sort(a, a + 3, [] (CanvasPoint a, CanvasPoint b) {
		return a.y < b.y;
		});
	CanvasPoint top = a[0];
	CanvasPoint mid = a[1];
	CanvasPoint bottom = a[2];
	
	// divide to 2 triangles (find 4th ver)
	CanvasPoint right, left, n = fourthPoint(top, mid, bottom);
	// find right and left points
	if (n.x > mid.x) {
		right = n;
		left = mid;
	}
	else {
		right = mid;
		left = n;
	}
	
	// color top triangle
	std::vector<float> rightSide = interpolateSingleFloats(top.x, right.x, int(right.y) - int(top.y));
	std::vector<float> leftSide = interpolateSingleFloats(top.x, left.x, int(left.y) - int(top.y));
	for (int i = top.y; i < int(right.y); i++) {
		CanvasPoint from = CanvasPoint(leftSide.at(i - top.y), i);
		CanvasPoint to = CanvasPoint(rightSide.at(i - top.y), i);
		drawLine(from, to, window, colour);
	}
	// color bottom triangle
	rightSide = interpolateSingleFloats(right.x, bottom.x, int(bottom.y) - int(right.y));
	leftSide = interpolateSingleFloats(left.x, bottom.x, int(bottom.y) - int(left.y));
	for (int i = mid.y; i < int(bottom.y); i++) {
		CanvasPoint from = CanvasPoint(leftSide.at(i - mid.y), i);
		CanvasPoint to = CanvasPoint(rightSide.at(i - mid.y), i);
		drawLine(from, to, window, colour);
	}
}

void drawPtCloud(DrawingWindow& window) {
	std::vector<ModelTriangle> modelTriangles = parseObj(OBJ_PATH);
	uint32_t colour = (255 << 24) + (255 << 16) + (255 << 8) + 255;
	for (auto mt : modelTriangles) {
		for (int i = 0; i < 3; i++) {
			CanvasPoint point = projectVertexOntoCanvasPoint(glm::vec3(0.0, 0.0, 4.0), mt.vertices[i], 2, camMat);
			window.setPixelColour(round(point.x), round(point.y), colour);
			std :: cout << "done" << std :: endl;
		}
	}
}

void drawWF(DrawingWindow& window) {
	std::vector<ModelTriangle> modelTriangles = parseObj(OBJ_PATH);
	Colour colour = Colour(255, 255, 255);
	for (auto mt : modelTriangles) {
		Colour colour = Colour(255, 255, 255);
		CanvasPoint v1 = projectVertexOntoCanvasPoint(glm::vec3(0.0, 0.0, 4.0), mt.vertices[0], 2, camMat);
		CanvasPoint v2 = projectVertexOntoCanvasPoint(glm::vec3(0.0, 0.0, 4.0), mt.vertices[1], 2, camMat);
		CanvasPoint v3 = projectVertexOntoCanvasPoint(glm::vec3(0.0, 0.0, 4.0), mt.vertices[2], 2, camMat);
		CanvasTriangle t = CanvasTriangle(v1, v2, v3);
		drawTriangle(t, window, colour);
	}
}

void drawObj(DrawingWindow& window, glm::vec3 camPos, glm::mat3 camMat) {
	std::vector<ModelTriangle> modelTriangles = parseObj(OBJ_PATH);
	Colour colour;
	for (auto mt : modelTriangles) {
		Colour colour = mt.colour;
		CanvasPoint v1 = projectVertexOntoCanvasPoint(camPos, mt.vertices[0], 2, camMat);
		CanvasPoint v2 = projectVertexOntoCanvasPoint(camPos, mt.vertices[1], 2, camMat);
		CanvasPoint v3 = projectVertexOntoCanvasPoint(camPos, mt.vertices[2], 2, camMat);
		CanvasTriangle t = CanvasTriangle(v1, v2, v3);
		drawFilTri(t, window, colour);
	}
}

void drawColour(DrawingWindow &window) {
	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow
	window.clearPixels();

	std::vector<glm::vec3> valLeft;
	std::vector<glm::vec3> valRight;
	valLeft = interpolateThreeElementValues(topLeft, bottomLeft, window.width);
	valRight = interpolateThreeElementValues(topRight, bottomRight, window.width);
	for (size_t y = 0; y < window.height; y++) {
		std::vector<glm::vec3> valLR;
		valLR = interpolateThreeElementValues(valLeft[y], valRight[y], window.width);
		for (size_t x = 0; x < window.width; x++) {
			float red = valLR[x].x;
			float green = valLR[x].y;
			float blue = valLR[x].z;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}

}

void draw(DrawingWindow &window){
	if (orbit) {
		camMat = rotateYaxis(window, camMat);
		camMat = lookAt(camPos, camMat);
	}

	if (mode == 1) {
		drawWF(window);
	}
	else if (mode == 2) {
		drawObj(window, camPos, camMat);
	}
	else if (mode == 3) {
		drawRayTraceAngle(window, camPos, camMat);
	}
		
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u){
			drawTriangle(randTri(), window, randCol());
			}
		else if (event.key.keysym.sym == SDLK_f){
			drawFilTri(randTri(), window, randCol());
		}
		else if (event.key.keysym.sym == SDLK_a) {
			resetWindow(window);
			camPos.x += 0.1;}
		else if (event.key.keysym.sym == SDLK_d) {
			resetWindow(window);
			camPos.x -= 0.1;}
		else if (event.key.keysym.sym == SDLK_w){camPos = moveUp(window, camPos);}
		else if (event.key.keysym.sym == SDLK_s){camPos = moveDown(window, camPos);}
		else if (event.key.keysym.sym == SDLK_e) {
			resetWindow(window);
			camPos.z -= 0.1;}
		else if (event.key.keysym.sym == SDLK_q) {
			resetWindow(window);
			camPos.z += 0.1;}
		else if (event.key.keysym.sym == SDLK_x){camMat = rotateXaxis(window, camMat);}
		else if (event.key.keysym.sym == SDLK_y){camMat = rotateYaxis(window, camMat);}
		else if (event.key.keysym.sym == SDLK_o) orbit = orbit ? false : true;

		else if (event.key.keysym.sym == SDLK_1) mode = 1;
		else if (event.key.keysym.sym == SDLK_2) mode = 2;
		else if (event.key.keysym.sym == SDLK_3) mode = 3;
		
		
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
