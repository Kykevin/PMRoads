#pragma once


#include <iostream>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "GL/glut.h"


#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "glm/vec4.hpp"
#include "glm/mat4x4.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

using std::cout;
using std::endl;
using std::vector;

using namespace glm;

struct Road{
	vec2 startVertex;
	vec2 endVertex;
	vector<vec2> middleVertex;
	vec2 centerVertex;
	float lineWidth;
};

struct Point{
	vec2 pos;
	int count;
};

/* Constant values */
const float PI = 3.14159265359;

/* Window information */
float windowWidth = 1200;
float windowHeight = 900;

float imageWidth = 600;
float imageHeight = 600;

/* #### The modelview matrix to load; remember OpenGL is row-major and not column major */
/* You do not necessarily need to use this matrix */
GLfloat myModelMat[4][4] = {
	{ 1, 0, 0, 0 },
	{ 0, 1, 0, 0 },
	{ 0, 0, 1, 0 },
	{ -5, -20, -26, 1 }
};

/* #### Define variables for the camera here... */

/* Field of view, aspect and near/far planes for the perspective projection */
float fovy = 45.0;
float aspect = windowWidth / windowHeight;
float zNear = 0.1;
float zFar = 1000.0;

float tx = 0;
float ty = 0;
float tz = -3;
float rx = 62;
float ry = 0;
float rz = 0;
float rotateFactor = 2;

bool isPerspective = true;

/* Vertices for the floor plane */
float floorVertices[12] = 
{
	0.0,  0.0, 0.0,
	0.0,  0.0,  10.0,
	10.0,  0.0,  10.0,
	10.0,  0.0, 0.0
};



vector<vec2> vertices;
vector<Road> roads;
vector<vector<Road>> largeAreas;
vector<vector<Road>> tempAreas;
vector<float> areaAngle;
vector<vec2> tempCenters;

float currxmax,currxmin,currymax,currymin,currArea;

float lineColors[6] = {
	1.0, 228/255.0, 181.0/255.0,
	1.0, 228/255.0, 181.0/255.0
};

/* IDs of the menu and the submenus */
int menuID;

//image data
//unsigned char imageData[600][600];
unsigned char reds[600][600];
unsigned char greens[600][600];
unsigned char blues[600][600];
int rows;
int cols;

int depth = 0;



float mapToPlane(float x){
	return x * 10;
}

vec2 findPerpNorm(const vec2& vec){
	return vec2(-normalize(vec).y,normalize(vec).x);
}


int orientation(vec2 p0,vec2 p1,vec2 p2)
{
	//calculate cross product of (p1-p0)X (p2-p0)

	int x1= p1.x- p0.x;
	int y1 = p1.y - p0.y;
	int x2= p2.x- p0.x;
	int y2 = p2.y - p0.y;

	int determinant = x1*y2 - x2*y1;

	if(determinant==0)
		return 0;

	if(determinant>0)
		return 1;

	if(determinant<0)
		return 2;

}

int onsegment(vec2 p0,vec2 p1,vec2 p2)
{
	if(p2.x >= min(p0.x,p1.x)&& p2.x <= max(p0.x,p1.x))
		return 1;

	return 0;
}

int dointersect(vec2 p1,vec2 q1,vec2 p2,vec2 q2)
{
	int o1 = orientation(p1,q1,q2);
	int o2 = orientation(p1,q1,p2);
	int o3 = orientation(p2,q2,p1);
	int o4 = orientation(p2,q2,q1);

	if(o1!=o2&&o3!=o4)  //handles general cases
		return 1;

	if(o1==0&&o2==0&&o3==0&&o4==0)  //handles special cases when all four points are collinear
	{
		if(onsegment(p1,q1,p2)||onsegment(p1,q1,q2))
			return 1;

	}
	return 0;

}

vec2 findIntersectPoint(vec2 start1, vec2 end1, vec2 start2, vec2 end2){
	float k1, k2;
	float b1, b2;
	if(abs(start1.x - end1.x) <= 0.00001){
		k2 = (start2.y-end2.y)/(start2.x-end2.x);
		return vec2(start1.x,k2*(start1.x-start2.x)+start2.y);
	}
	if(abs(start2.x - end2.x) <= 0.00001){
		k1 = (start1.y-end1.y)/(start1.x-end1.x);
		return vec2(start2.x,k1*(start2.x-start1.x)+start1.y);
	}
	k1 = (start1.y-end1.y)/(start1.x-end1.x);
	k2 = (start2.y-end2.y)/(start2.x-end2.x);
	b1 = start1.y -	k1 * start1.x;
	b2 = start2.y -	k2 * start2.x;
	//cout << "find intersect" << endl;
	//cout << start1.x << " " << start1.y << endl;
	//cout << end1.x << " " << end1.y << endl;
	//cout << start2.x << " " << start2.y << endl;
	//cout << end2.x << " " << end2.y << endl;
	//cout << "end" << endl;
	return vec2((b2-b1)/(k1-k2), (b2-b1)/(k1-k2)*k1+b1);
}

float crossProduct(vec2 line1, vec2 line2){
	return line1.x*line2.y - line1.y*line2.x;
}

float findDensity(vector<Road>& currRoads, float currAngle){
	float xmin = 9999,xmax = -999,ymin = 9999,ymax = -999;
	bool first = true;
	mat2 rotation = mat2(cos(currAngle),-sin(currAngle),sin(currAngle),cos(currAngle));
	currArea = 999999;
	for(vector<Road>::iterator it = currRoads.begin();it != currRoads.end();++it){
		vec2 temp = rotation * (*it).startVertex;
		xmin = xmin > temp.x ? temp.x : xmin;
		xmax = xmax < temp.x ? temp.x : xmax;
		ymin = ymin > temp.y ? temp.y : ymin;
		ymax = ymax < temp.y ? temp.y : ymax;
		temp = rotation * (*it).endVertex;
		xmin = xmin > temp.x ? temp.x : xmin;
		xmax = xmax < temp.x ? temp.x : xmax;
		ymin = ymin > temp.y ? temp.y : ymin;
		ymax = ymax < temp.y ? temp.y : ymax;
	}
	currArea = (xmin-xmax)*(ymin-ymax);
	currxmax = xmax;
	currxmin = xmin;
	currymin = ymin;
	currymax = ymax;

	float density = 0;
	float count = 0;
	rotation = mat2(cos(-currAngle),-sin(-currAngle),sin(-currAngle),cos(-currAngle));
	for(float i = currxmin; i <= currxmax; i+= (currxmax-currxmin)/10){
		for(float j = currymin; j <= currymax; j+= (currymax-currymin)/10){
			vec2 temp = rotation * vec2(i,j);
			density += reds[(int)temp.y][(int)temp.x];
			count++;
		}
	}
	//cout << "density: " << density/count << endl;
	return density/count;
}

vector<Road> sortLoop(vector<Road>& loop){
	vector<Road> newLoop;
	newLoop.push_back(loop.back());
	loop.pop_back();
	vec2 lastEnd = newLoop.front().endVertex;
	while(!loop.empty()){
		for(vector<Road>::iterator it = loop.begin();it != loop.end();++it){
			if((*it).startVertex == lastEnd || (*it).endVertex == lastEnd){
				if((*it).endVertex == lastEnd){
					vec2 temp = (*it).startVertex;
					(*it).startVertex = (*it).endVertex;
					(*it).endVertex = temp;
				}
				lastEnd = (*it).endVertex;
				newLoop.push_back(*it);
				loop.erase(it);
				break;
			}
		}
	}
	return newLoop;
}


void generateOuterLoop(){
	int i,j;
	int i2,j2;
	vec2 start(0,0);
	for(i = 5; i < 595; ++i){
		if(start.x != 0) break;
		for(j = 5; j < 595;++j){
			if(start.x != 0) break;
			bool chk = true;
			for(i2 = -3; i2 < 4; i2++){
				for(j2 = -3; j2 < 4; j2++){
					if(blues[i+i2][j+j2] >= 20){
						chk = false;
					}
				}
			}
			if(chk){
				start.y = i;
				start.x = j;
			}
		}
	}
	vertices.push_back(start);
	vec2 currPoint = start;

	float angleStep = 3/180.0*PI;
	float angle = 0;
	float lastAngle = 0;
	mat2 rotation(cos(angle),-sin(angle),sin(angle),cos(angle));
	mat2 lastRotation;
	vec2 step1(50,0),step2(100,0),step3(150,0);
	bool chk1,chk2,chk3;
	bool lastChk1,lastChk2,lastChk3;
	bool first = true;
	int result;
	//rotate left
	int counter = 0;
	while(1){
		if(distance(start,currPoint) < 80 && counter >= 6){
			Road road;
			road.centerVertex = vec2(-1,-1);
			road.endVertex = start;
			road.startVertex = currPoint;
			road.lineWidth = 2;
			roads.push_back(road);
			break;
		}
		if(first != true){
			lastChk1 = chk1;
			lastChk2 = chk2;
			lastChk3 = chk3;
			chk1 = false;
			chk2 = false;
			chk3 = false;
			lastRotation = rotation;
			rotation = mat2(cos(angle),-sin(angle),sin(angle),cos(angle));
			vec2 temp = currPoint + rotation * step1;
			if((int)temp.y >= 0 && (int)temp.y < 600 && (int)temp.x >= 0 && (int)temp.x < 600){
				if(blues[(int)temp.y][(int)temp.x] <= 20) {
					chk1 = true;
				}
			}
			temp = currPoint + rotation * step2;
			if((int)temp.y >= 0 && (int)temp.y < 600 && (int)temp.x >= 0 && (int)temp.x < 600){
				if(blues[(int)temp.y][(int)temp.x] <= 20) {
					chk2 = true;
				}
			}
			temp = currPoint + rotation * step3;
			if((int)temp.y >= 0 && (int)temp.y < 600 && (int)temp.x >= 0 && (int)temp.x < 600){
				if(blues[(int)temp.y][(int)temp.x] <= 20) {
					chk3 = true;
				}
			}
			//cout << (int)temp.y << "   " << (int)temp.x << endl;

			if(chk1 == false && lastChk1 == true){
				counter ++;
				if(chk3 == false && lastChk3 == true){
					temp = currPoint + lastRotation * step3;
					lastAngle = angle-angleStep;
					vertices.push_back(temp);
					Road road;
					road.centerVertex = vec2(-1,-1);
					road.endVertex = temp;
					road.startVertex = currPoint;
					road.lineWidth = 2;
					roads.push_back(road);
					currPoint = temp;
				}
				else if(chk2 == false && lastChk2 == true){
					temp = currPoint + lastRotation * step2;
					lastAngle = angle-angleStep;
					vertices.push_back(temp);
					Road road;
					road.centerVertex = vec2(-1,-1);
					road.endVertex = temp;
					road.startVertex = currPoint;
					road.lineWidth = 2;
					roads.push_back(road);
					currPoint = temp;
				}
				else {
					temp = currPoint + lastRotation * step1;
					lastAngle = angle-angleStep;
					vertices.push_back(temp);
					Road road;
					road.centerVertex = vec2(-1,-1);
					road.endVertex = temp;
					road.startVertex = currPoint;
					road.lineWidth = 2;
					roads.push_back(road);
					currPoint = temp;
				}
				angle = angle - 20.0 * angleStep;
			}

		} 
		else {
			rotation = mat2(cos(angle),-sin(angle),sin(angle),cos(angle));
			vec2 temp = currPoint + rotation * step1;
			chk1 = false;
			chk2 = false;
			chk3 = false;
			if((int)temp.y >= 0 && (int)temp.y < 600 && (int)temp.x >= 0 && (int)temp.x < 600){
				if(blues[(int)temp.y][(int)temp.x] <= 20) {
					chk1 = true;
				}
			}
			temp = currPoint + rotation * step2;
			if((int)temp.y >= 0 && (int)temp.y < 600 && (int)temp.x >= 0 && (int)temp.x < 600){

				if(blues[(int)temp.y][(int)temp.x] <= 20) {
					chk2 = true;
				}
			}
			temp = currPoint + rotation * step3;
			if((int)temp.y >= 0 && (int)temp.y < 600 && (int)temp.x >= 0 && (int)temp.x < 600){

				if(blues[(int)temp.y][(int)temp.x] <= 20) {
					chk3 = true;
				}
			}
			first = false;
		}
		angle += angleStep;
	}
}


void generateLargeArea(vector<Road>& currRoads, int type){

	//depth++;

	//if( depth >= 2) {
	//	depth--;
	//	//return;
	//}

	float angle = 0;
	float angleStep = 3*PI/180.0;
	float minArea = 999999;
	vec2 centerVector;
	vec2 perpVector;
	vec2 centerPoint;
	mat2 rotation;
	float tempAngle;
	for(angle = 0;angle < 2*PI; angle+=angleStep){
		rotation = mat2(cos(angle),-sin(angle),sin(angle),cos(angle));
		float xmin = 9999,xmax = -999,ymin = 9999,ymax = -999;
		for(vector<Road>::iterator it = currRoads.begin();it != currRoads.end();++it){
			vec2 temp = rotation * (*it).startVertex;
			xmin = xmin > temp.x ? temp.x : xmin;
			xmax = xmax < temp.x ? temp.x : xmax;
			ymin = ymin > temp.y ? temp.y : ymin;
			ymax = ymax < temp.y ? temp.y : ymax;
			temp = rotation * (*it).endVertex;
			xmin = xmin > temp.x ? temp.x : xmin;
			xmax = xmax < temp.x ? temp.x : xmax;
			ymin = ymin > temp.y ? temp.y : ymin;
			ymax = ymax < temp.y ? temp.y : ymax;
		}
		if((xmin-xmax)*(ymin-ymax) < minArea){
			minArea = (xmin-xmax)*(ymin-ymax);
			centerPoint.x = (xmax + xmin)/2;
			centerPoint.y = (ymax + ymin)/2;
			centerPoint = mat2(cos(-angle),-sin(-angle),sin(-angle),cos(-angle)) * centerPoint;
			if( xmax - xmin > ymax - ymin){
				perpVector = rotation * vec2(0,1);
				centerVector = rotation * vec2(xmax-xmin,0);
				tempAngle = angle;
			}
			else{
				perpVector = rotation * vec2(1,0);
				centerVector = rotation * vec2(0,ymax - ymin);
				tempAngle = angle;
			}
		}
	}
	centerPoint = centerPoint + ((rand() % 10 - 5) / 50.0f)*centerVector;
	if(minArea <= 25000 && type == 0) {
		//depth--;
		//cout << "minArea: " << minArea << "   Angle: " << tempAngle<< endl;
		currRoads = sortLoop(currRoads);
		largeAreas.push_back(currRoads);
		areaAngle.push_back(tempAngle);
		return;
	}
	if(minArea <= 80000 && type == 1) {
		//int pick = rand()%2;
		//if(pick == 0){
		currRoads = sortLoop(currRoads);
		tempAreas.push_back(currRoads);
		tempCenters.push_back(centerPoint);
		return;
		//}
	}

	//cout << "minArea: " << minArea << endl;
	//cout << "angle: " << tempAngle << endl;
	//cout << "centerPoint: " << centerPoint.x << "   " << centerPoint.y << endl;
	//cout << "centerVector: " << centerVector.x << "   " << centerVector.y << endl;
	//cout << "perpVec: " << perpVector.x << "   " << perpVector.y << endl;



	//angle = (rand()%10-5)*PI/180.0; 
	//perpVector = mat2(cos(angle),-sin(angle),sin(angle),cos(angle)) * perpVector;
	vector<vec2> newVecs;
	vector<Road> del,add;
	for(vector<Road>::iterator it = currRoads.begin();it != currRoads.end();++it){
		vec2 tempInter = findIntersectPoint(centerPoint,centerPoint+perpVector,(*it).startVertex,(*it).endVertex);
		//if(abs((*it).startVertex.y-(*it).endVertex.y) <= 0.0001){
		//	cout << "line " << (*it).startVertex.x << "   " <<(*it).startVertex.y << endl;		
		//	cout << "line " << (*it).endVertex.x << "   " <<(*it).endVertex.y << endl;		
		//	cout << "point " << tempInter.x << "  " << tempInter.y << endl;
		//}
		if((tempInter.x <= (*it).startVertex.x + 0.005 && tempInter.x >= (*it).endVertex.x - 0.005)
			|| (tempInter.x <= (*it).endVertex.x  + 0.005 && tempInter.x > (*it).startVertex.x- 0.005)){
				if((tempInter.y <= (*it).startVertex.y + 0.005 && tempInter.y >= (*it).endVertex.y- 0.005) 
					||(tempInter.y <= (*it).endVertex.y  + 0.005 && tempInter.y >= (*it).startVertex.y- 0.005)){
						//cout << tempInter.x << "  " << tempInter.y << endl;
						if(distance((*it).startVertex,(*it).endVertex) <= sqrt(minArea) / 4.0f){
							if(distance((*it).startVertex,tempInter) < distance((*it).startVertex,tempInter)){
								vec2 vertex1 = (*it).startVertex;
								newVecs.push_back(vertex1);
							}
							else{
								vec2 vertex1 = (*it).endVertex;
								newVecs.push_back(vertex1);
							}
						}				
						else {
							vec2 vertex1;
							vertex1 = tempInter;
							//vertices.push_back(vertex1);
							(*it).middleVertex.push_back(vertex1);

							newVecs.push_back(vertex1);
							Road newRoad1;
							newRoad1.centerVertex = vec2(-1,-1);
							newRoad1.startVertex = vertex1;
							newRoad1.endVertex = (*it).startVertex;
							newRoad1.lineWidth = 2;

							Road newRoad2;
							newRoad2.centerVertex = vec2(-1,-1);
							newRoad2.startVertex = vertex1;
							newRoad2.endVertex = (*it).endVertex;
							newRoad2.lineWidth = 2;

							add.push_back(newRoad1);
							add.push_back(newRoad2);
							del.push_back(*it);
						}

				}

		}
	}

	for(vector<Road>::iterator it = del.begin();it != del.end();++it){
		for(vector<Road>::iterator it2 = currRoads.begin();it2 != currRoads.end();++it2){
			if((*it).startVertex == (*it2).startVertex && (*it).endVertex == (*it2).endVertex){
				currRoads.erase(it2);
				break;
			}
		}
	}
	for(vector<Road>::iterator it = add.begin();it != add.end();++it){
		currRoads.push_back((*it));
	}

	if(newVecs.size() != 2){
		//largeAreas.push_back(currRoads);
		return;
	}

	vec2 vertex1 = newVecs.back();
	newVecs.pop_back();
	vec2 vertex2 = newVecs.back();
	newVecs.pop_back();
	//cout << "vertex1 and 2: ";
	//cout << vertex1.x <<"   "<< vertex1.y <<  "----" << vertex2.x <<"   "<< vertex2.y << endl ;
	Road newRoad;
	newRoad.centerVertex = vec2(-1,-1);
	newRoad.endVertex = vertex2;
	newRoad.startVertex = vertex1;
	newRoad.lineWidth = 2;
	roads.push_back(newRoad);
	currRoads.push_back(newRoad);
	vector<Road> nextRoads1;
	vector<Road> nextRoads2;
	for(vector<Road>::iterator it = currRoads.begin();it != currRoads.end();++it){
		if(crossProduct( (*it).startVertex-vertex1,vertex1-vertex2) < 0 || crossProduct( (*it).endVertex-vertex1,vertex1-vertex2) < 0){
			nextRoads1.push_back(*it);
		}
		if(crossProduct( (*it).startVertex-vertex1,vertex1-vertex2) > 0 || crossProduct( (*it).endVertex-vertex1,vertex1-vertex2) > 0){
			nextRoads2.push_back(*it);
		}
	}
	nextRoads1.push_back(newRoad);
	nextRoads2.push_back(newRoad);

	generateLargeArea(nextRoads1,type);
	generateLargeArea(nextRoads2,type);

	//depth--;
}

struct myLineLesser
{
	bool operator()( const vec2& lx, const vec2& rx ) const {
		return lx.x < rx.x;
	}
};

vector<vec2> currRoadIntersectOuterLoop(vector<Road>& currRoads, vec2 startPoint, vec2 endPoint){
	vector<vec2> newVecs;
	for(vector<Road>::iterator it = currRoads.begin();it != currRoads.end();++it){
		vec2 tempInter = findIntersectPoint(startPoint,endPoint,(*it).startVertex,(*it).endVertex);
		if((tempInter.x <= (*it).startVertex.x + 0.005 && tempInter.x >= (*it).endVertex.x - 0.005)
			|| (tempInter.x <= (*it).endVertex.x  + 0.005 && tempInter.x > (*it).startVertex.x- 0.005)){
				if((tempInter.y <= (*it).startVertex.y + 0.005 && tempInter.y >= (*it).endVertex.y- 0.005) 
					||(tempInter.y <= (*it).endVertex.y  + 0.005 && tempInter.y >= (*it).startVertex.y- 0.005)){
						vec2 vertex1;
						vertex1 = tempInter;
						newVecs.push_back(vertex1);
				}
		}
	}
	return newVecs;

}

void generateGrid(vector<Road>& currRoads, float currAngle){
	float density = findDensity(currRoads, currAngle);
	float shortEdgeStep = 15 - density * 7.0f/255.0f;
	float longEdgeStep = 30 - density * 15.0f/255.0f;
	float xstep, ystep;

	if(currxmax - currxmin < currymax - currymin){
		shortEdgeStep = (int)((currxmax-currxmin) / shortEdgeStep);
		longEdgeStep = (int)((currymax-currymin) / longEdgeStep);
		xstep = (currxmax-currxmin) / shortEdgeStep;
		ystep = (currymax-currymin) / longEdgeStep;
	}
	else{
		shortEdgeStep = (int)((currymax-currymin) / shortEdgeStep);
		longEdgeStep = (int)((currxmax-currxmin) / longEdgeStep);
		ystep = (currymax-currymin) / shortEdgeStep;
		xstep = (currxmax-currxmin) / longEdgeStep;
	}
	//cout << "xstep: " << xstep << "    ystep: " << ystep << endl;
	//mat2 rotation(cos(-currAngle),-sin(-currAngle),sin(-currAngle),cos(-currAngle));
	////cout << "currArea: " << currArea << "   currAngle:" << currAngle <<endl;
	//vec2 a = rotation * vec2(currxmin,currymin);
	//vec2 b = rotation * vec2(currxmax,currymax);
	//cout << "currxmin: " << a.x << "   ";
	//cout << "currxmax: " << b.x << "   ";
	//cout << "currymin: " << a.y << "   ";
	//cout << "currymax: " << b.y << "   ";
	//cout << endl;
	mat2 rotation = mat2(cos(-currAngle),-sin(-currAngle),sin(-currAngle),cos(-currAngle));
	vec2 xvector = rotation*vec2(1,0);
	vec2 yvector = rotation*vec2(0,1);


	for(float currx = currxmin + xstep; currx < currxmax; currx+=xstep){
		vec2 temp = rotation*vec2(currx,currymin);
		//cout << "linea " << temp.x << "   " << temp.y << endl;		
		//cout << "lineb " << (temp+yvector).x << "   " << (temp+yvector).y << endl;		
		vector<vec2> newVecs = currRoadIntersectOuterLoop(currRoads,temp, temp+yvector);
		if(newVecs.size()%2 == 0){
			sort(newVecs.begin(),newVecs.end(),myLineLesser());
			while(!newVecs.empty()){
				vec2 vertex1 = newVecs.back();
				newVecs.pop_back();
				vec2 vertex2 = newVecs.back();
				newVecs.pop_back();
				//cout << "vertex1 and 2: ";
				//cout << vertex1.x <<"   "<< vertex1.y <<  "----" << vertex2.x <<"   "<< vertex2.y << endl ;
				Road newRoad;
				newRoad.centerVertex = vec2(-1,-1);
				newRoad.endVertex = vertex2;
				newRoad.startVertex = vertex1;
				newRoad.lineWidth = 1;
				roads.push_back(newRoad);
				//cout << "succeed once" << endl;
			}
		}
	}
	for(float curry = currymin + ystep; curry < currymax; curry+=ystep){
		vec2 temp = rotation*vec2(currxmin,curry);
		//cout << "line " << temp.x << "   " << temp.y << endl;		
		//cout << "line " << (temp+xvector).x << "   " << (temp+xvector).y << endl;		
		vector<vec2> newVecs = currRoadIntersectOuterLoop(currRoads,temp, temp+xvector);
		if(newVecs.size()%2 == 0){
			sort(newVecs.begin(),newVecs.end(),myLineLesser());
			while(!newVecs.empty()){
				vec2 vertex1 = newVecs.back();
				newVecs.pop_back();
				vec2 vertex2 = newVecs.back();
				newVecs.pop_back();
				//cout << "vertex1 and 2: ";
				//cout << vertex1.x <<"   "<< vertex1.y <<  "----" << vertex2.x <<"   "<< vertex2.y << endl ;
				Road newRoad;
				newRoad.centerVertex = vec2(-1,-1);
				newRoad.endVertex = vertex2;
				newRoad.startVertex = vertex1;
				newRoad.lineWidth = 1;
				roads.push_back(newRoad);
				//cout << "succeed once" << endl;
			}
		}
	}
}

void generateRadial(vector<Road>& currRoads, vec2 centerPoint){
	vector<vec2> rays;
	vector<Road> rayLines;
	int pick = rand()%10;
	int degree = rand()%360;
	float angle;
	if(pick <= 6){
		for(int i = 0; i < 5;++i){
			degree = (degree + 60 + rand()%12) % 360;
			angle = degree * PI / 180.0f;
			rays.push_back(mat2(cos(angle),-sin(angle),sin(angle),cos(angle)) * vec2(1,0));
		}
	} 
	else if(pick <= 8){
		for(int i = 0; i < 5;++i){
			degree = (degree + 49 + rand()%12) % 360;
			angle = degree * PI / 180.0f;
			rays.push_back(mat2(cos(angle),-sin(angle),sin(angle),cos(angle)) * vec2(1,0));
		}
	}
	else{
		for(int i = 0; i < 5;++i){
			degree = (degree + 38 + rand()%12) % 360;
			angle = degree * PI / 180.0f;
			rays.push_back(mat2(cos(angle),-sin(angle),sin(angle),cos(angle)) * vec2(1,0));
		}
	}
	vector<Road> add;
	vector<Road> del;
	for(vector<vec2>::iterator rayIt = rays.begin();rayIt != rays.end();++rayIt){
		vector<vec2> newVecs;
		vec2 startPoint = centerPoint;
		vec2 endPoint = centerPoint+(*rayIt);
		vector<Road> interRoads;
		for(vector<Road>::iterator it = currRoads.begin();it != currRoads.end();++it){
			vec2 tempInter = findIntersectPoint(startPoint,endPoint,(*it).startVertex,(*it).endVertex);
			if((tempInter.x <= (*it).startVertex.x + 0.005 && tempInter.x >= (*it).endVertex.x - 0.005)
				|| (tempInter.x <= (*it).endVertex.x  + 0.005 && tempInter.x > (*it).startVertex.x- 0.005)){
					if((tempInter.y <= (*it).startVertex.y + 0.005 && tempInter.y >= (*it).endVertex.y- 0.005) 
						||(tempInter.y <= (*it).endVertex.y  + 0.005 && tempInter.y >= (*it).startVertex.y- 0.005)){
							vec2 vertex1;
							vertex1 = tempInter;
							newVecs.push_back(vertex1);
							interRoads.push_back(*it);
					}
			}
		}
		vec2 vertex;
		Road interRoad;
		while(!newVecs.empty()){
			vertex = newVecs.back();
			newVecs.pop_back();
			interRoad = interRoads.back();
			interRoads.pop_back();
			if((vertex.x - centerPoint.x) / (*rayIt).x >= 0){
				break;
			}
		}

		Road newRoad;
		newRoad.centerVertex = vec2(-1,-1);
		newRoad.endVertex = vertex;
		newRoad.startVertex = centerPoint;
		newRoad.lineWidth = 2;
		roads.push_back(newRoad);
		rayLines.push_back(newRoad);

		for(vector<Road>::iterator it2 = currRoads.begin();it2 != currRoads.end();++it2){
			if(interRoad.startVertex == (*it2).startVertex && interRoad.endVertex == (*it2).endVertex){
				Road newRoad1;
				newRoad1.centerVertex = vec2(-1,-1);
				newRoad1.startVertex = vertex;
				newRoad1.endVertex = interRoad.startVertex;
				newRoad1.lineWidth = 2;

				Road newRoad2;
				newRoad2.centerVertex = vec2(-1,-1);
				newRoad2.startVertex = vertex;
				newRoad2.endVertex = interRoad.endVertex;
				newRoad2.lineWidth = 2;
				currRoads.erase(it2);
				currRoads.push_back(newRoad1);
				currRoads.push_back(newRoad2);
				break;
			}
		}
	}

	currRoads = sortLoop(currRoads);

	vector<vector<Road>> newLoops;
	bool newStart = false;
	vector<Road>::iterator startPos = currRoads.end();
	vector<Road>* currNewLoop = NULL;
	vector<Road>::iterator it = currRoads.begin();
	//cout << "Total: " << currRoads.size() << endl;
	while(it != startPos){
		//if(newStart){
		//	newStart = false;
		//}
		for(vector<Road>::iterator rayIt = rayLines.begin();rayIt != rayLines.end();++rayIt){
			if((*it).startVertex == (*rayIt).endVertex){
				if(startPos == currRoads.end()){
					startPos = it;
				}
				if(currNewLoop != NULL){
					//cout << "curr: " << currNewLoop->size() << endl;
					newLoops.push_back(*currNewLoop);
				}
				currNewLoop = new vector<Road>;
				break;
			}
		}
		if(startPos != currRoads.end()){
			if(currNewLoop != NULL){
				currNewLoop->push_back(*it);
			}
		}
		++it;
		if(it == currRoads.end()){
			it = currRoads.begin();
		}
		if(it == startPos){
			//	cout << "curr: " << currNewLoop->size() << endl;
			newLoops.push_back(*currNewLoop);
		}
	}
	while(!newLoops.empty()){
		vector<Road> currLoop = newLoops.back();
		newLoops.pop_back();
		//cout << "curr: " << currLoop.size() << endl;
		vector<Road>::iterator ray1,ray2;
		for(vector<Road>::iterator rayIt = rayLines.begin();rayIt != rayLines.end();++rayIt){
			if(currLoop.front().startVertex == (*rayIt).endVertex){
				ray1 = rayIt;
			}
			if(currLoop.back().endVertex == (*rayIt).endVertex){
				ray2 = rayIt;
			}
		}
		currLoop.push_back(*ray1);
		currLoop.push_back(*ray2);
		//cout << "curr: " << currLoop.size() << endl;
		generateLargeArea(currLoop, 0);
	}
	//cout << "CenterX " << centerPoint.x << "    CenterY " << centerPoint.y << endl;
}


/**
*    Function invoked for drawing using OpenGL
*/
void display()
{
	static int frameCount=0;

	/* #### frame count, might come in handy for animations */
	if(frameCount == 30){
		generateOuterLoop();
	}
	else if(frameCount == 60){
		vector<Road> currRoads;
		vector<vec2> currVertices;
		for(vector<Road>::iterator it = roads.begin(); it != roads.end();++it){
			currRoads.push_back(*it);
		}
		for(vector<vec2>::iterator it = vertices.begin(); it != vertices.end();++it){
			currVertices.push_back(*it);
		}
		generateLargeArea(currRoads,1);
		//		cout << "Area count: " << largeAreas.size() << endl;
	}
	else if(frameCount >= 100 && frameCount % 30 == 0 && !tempAreas.empty()){
		//cout << "left: " << largeAreas.size() << endl;
		if(rand() % 5 >= 2){
			generateLargeArea(tempAreas.back(),0);
		}
		else {
			generateRadial(tempAreas.back(), tempCenters.back());
		}
		tempAreas.pop_back();
		tempCenters.pop_back();
		//frameCount = 50001;
		if(tempAreas.empty()){
			frameCount = 9900;
		}
	}
	else if(frameCount >= 10000 && frameCount % 15 == 0){
		//cout << "left: " << largeAreas.size() << endl;
		if(!largeAreas.empty()){
			generateGrid(largeAreas.back(), areaAngle.back());
			largeAreas.pop_back();
			areaAngle.pop_back();
		}
	}

	if(frameCount <= 50000){
		frameCount++;
	}

	/* Clear the window */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	/* Set the perspective projection */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (isPerspective) {
		gluPerspective(fovy, aspect, zNear, zFar);
	} 
	else {
		glOrtho(-8,8,-2,6,0.1,1000);
	}
	/* #### Load/set the model view matrix */
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glLoadMatrixf((GLfloat *)myModelMat);

	/* Enable client */
	glEnableClientState(GL_VERTEX_ARRAY);

	/* Draw the floor */
	glVertexPointer(3, GL_FLOAT, 0, floorVertices);
	glTranslatef(5, 0, 26);
	glRotatef(rx,1,0,0);
	glRotatef(ry,0,1,0);
	glTranslatef(-5, 0, -26);
	glTranslatef(tx, 0, tz);
	glColor3f(132.0 / 255, 203.0 / 255, 107.0 / 255);
	glDrawArrays(GL_POLYGON, 0, 4);

	vector<Road>::iterator it = roads.begin();
	for(it;it != roads.end();it++){
		float lineVertices[12] = {((*it).startVertex.x-(*it).lineWidth * findPerpNorm((*it).startVertex-(*it).endVertex).x)/60.0,
			0.01, ((*it).startVertex.y-(*it).lineWidth * findPerpNorm((*it).startVertex-(*it).endVertex).y)/60.0, 
			((*it).endVertex.x-(*it).lineWidth * findPerpNorm((*it).startVertex-(*it).endVertex).x)/60.0, 
			0.01, ((*it).endVertex.y-(*it).lineWidth * findPerpNorm((*it).startVertex-(*it).endVertex).y)/60.0,
			((*it).endVertex.x+(*it).lineWidth * findPerpNorm((*it).startVertex-(*it).endVertex).x)/60.0, 
			0.01, ((*it).endVertex.y+(*it).lineWidth * findPerpNorm((*it).startVertex-(*it).endVertex).y)/60.0,
			((*it).startVertex.x+(*it).lineWidth * findPerpNorm((*it).startVertex-(*it).endVertex).x)/60.0, 
			0.01, ((*it).startVertex.y+(*it).lineWidth * findPerpNorm((*it).startVertex-(*it).endVertex).y)/60.0};
		glVertexPointer(3, GL_FLOAT, 0, lineVertices);
		if((*it).lineWidth == 2){
			glColor3f(255.0/255.0, 228/255.0, 181.0/255.0);
		}
		else{
			glColor3f(255.0/255.0, 255/255.0, 224.0/255.0);
		}
		glDrawArrays(GL_POLYGON, 0, 4);
	}

	/* Disable client */
	glDisableClientState(GL_VERTEX_ARRAY);

	/* Force execution of OpenGL commands */
	glFlush();

	/* Swap buffers for animation */
	glutSwapBuffers();
}


/**
*    Function invoked when window system events are not being received
*/
void idle()
{
	/* Redraw the window */
	glutPostRedisplay();
}


/**
*    #### Function invoked when an event on regular keys occur
*/
void keyboard(unsigned char k, int x, int y)
{
	/* Show which key was pressed */
	//cout << "Pressed \"" << k << "\" ASCII: " << (int)k << endl;

	if(k == 'a') 
	{
		ry -= rotateFactor;
	}
	else if(k == 'd') 
	{
		ry += rotateFactor;
	}
	else if(k == 'w') 
	{
		rx -= rotateFactor;
	}
	else if(k == 's') 
	{
		rx += rotateFactor;
	}
	else if (k == 'o'){
		isPerspective = false;
	}
	else if (k == 'p'){
		isPerspective = true;
	}
	else if (k == 'r'){
		rx = 0;
		ry = 0;
		rz = 0;
		tx = 0;
		ty = 0;
		tz = 0;
	}
	else if (k == 27)
	{

		exit(0);
	}
}


/**
*	#### Function invoked when an event on special keys occur
*/
void special(int key, int x, int y) 
{
	if(key == GLUT_KEY_UP) 
	{
		tz += 0.5;
	}
	else if(key == GLUT_KEY_DOWN) 
	{
		tz -= 0.5;
	}
	else if(key == GLUT_KEY_RIGHT) 
	{		
		tx -= 0.5;
	}
	else if(key == GLUT_KEY_LEFT) 
	{
		tx += 0.5;
	}
	//printf("tx = %f\tty = %f\ttz = %f\n",tx,ty,tz);
}



/**
*	Function called when the an option of the menu is selected
*/
void menu(int value)
{
	if (value == 0)
	{
		exit(0);
	}
}


void makeMenu()
{

	menuID = glutCreateMenu(menu);

	glutAddMenuEntry("Exit", 0);

	/* Attach menu to the right click */
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void readImage(){

}

/**
*	Set OpenGL initial state
*/
void init()
{
	/* Set clear color */
	glClearColor(1.0, 1.0, 1.0, 0.0);
	//readImage();
}

void setTestCase(int num);

/**
*	Main function
*/

int main(int argc, char **argv)
{
	//srand(time(NULL));
	int num;
	cout << "Input test case number: ";
	std::cin >> num;

	setTestCase(num);

	/* Initialize the GLUT window */
	glutInit(&argc, argv);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(30, 30);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Final Project");


	/* Set OpenGL initial state */
	init();

	/* Callback functions */
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(special);
	makeMenu();

	/* Start the main GLUT loop */
	/* NOTE: No code runs after this */
	glutMainLoop();

	return 0;
}


void drawRect(int xmin, int xmax, int ymin, int ymax, char r, char b){
	int i,j;
	for(i = 0;i < 600;++i){
		for(j = 0;j < 600;++j){
			if(i >= xmin && i <= xmax){
				if(j >= ymin && j <= ymax){
					if(b != 0)
						blues[j][i] = b;
					if(r != 0)
						reds[j][i] = r;
				}
			}
		}
	}
}
void drawCircle(int x, int y, int rad, char r, char b){
	int i,j;
	for(i = 0;i < 600;++i){
		for(j = 0;j < 600;++j){
			if (pow(i-x,2) + pow(j-y,2) <= pow(rad,2)){
				if(b != 0)
					blues[j][i] = b;
				if(r != 0)
					reds[j][i] = r;
			}
		}
	}
}

void setTestCase(int num){
	int i,j;
	if(num == 1){
		for(i = 0;i < 600;++i){
			for(j = 0;j < 600;++j){
				blues[i][j] = 255;
				greens[i][j] = 0;
				if(i >= 30 && i <= 570){
					if(j >= 100 && j <= 500){
						greens[i][j] = 100;
						blues[i][j] = 0;
					}
				}
			}
		}

		for(i = 0;i < 600;++i){
			for(j = 0;j < 600;++j){
				if (pow(i-300,2) + pow(j-300,2) <= pow(250,2)){
					//blues[i][j] = 0;
					//greens[i][j] = 100;
				}
				else {
					greens[i][j] = 0;
					blues[i][j] = 255;
				}
			}
		}
	}
	if(num == 2){
		for(i = 0;i < 600;++i){
			for(j = 0;j < 600;++j){
				blues[i][j] = 255;
				greens[i][j] = 0;
				if(i >= 100 && i <= 500){
					if(j >= 30 && j <= 570){
						greens[i][j] = 100;
						blues[i][j] = 0;
					}
				}
			}
		}
		for(i = 0;i < 600;++i){
			for(j = 0;j < 600;++j){
				if(i >= 100 && i <= 300){
					if(j >= 100 && j <= 300){
						greens[i][j] = 0;
						reds[i][j] = 255;
					}
					if(j >= 300 && j <= 500){
						greens[i][j] = 0;
						reds[i][j] = 255;
					}
				}
			}
		}
	}
	if(num == 3){
		drawRect(0,600,0,600,0,255);
		drawRect(100,500,100,250,0,1);
		drawRect(150,450,250,350,0,1);
		drawCircle(300,425,125,0,1);
		drawCircle(300,450,100,255,0);
	}
	if(num == 4){
		drawRect(0,600,0,600,0,255);

		drawRect(100,500,100,250,0,1);
		drawRect(150,450,250,350,0,1);
		drawCircle(300,325,125,0,1);
		drawRect(100,500,400,550,0,1);


		drawCircle(300,350,100,255,0);
	}

	if(num == 5){
		drawRect(0,600,0,600,0,255);
		
		drawCircle(200,200,150,0,1);
		drawCircle(400,400,150,0,1);
		drawCircle(200,400,150,0,1);
		drawCircle(400,200,150,0,1);
		drawCircle(254,300,200,0,1);

		drawRect(200,400,200,400,125,0);
		drawCircle(300,350,100,255,0);	
	}
}