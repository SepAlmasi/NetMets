#include "rtsFiberNetwork.h"
#include "rts_glShaderProgram.h"
#include "GL/glut.h"
#include "rts_glutRenderWindow.h"
#include <time.h>

extern void ComputeNetMets();
extern rtsFiberNetwork* goldNetwork;
extern rtsFiberNetwork* testNetwork;
extern float sigmaG, sigmaC;
float network_span;
CoreGraphList coreGraph;
vector<point3D<float> > sequenceColors; 
int current_sequence = 0;

//shader variables
rts_glShaderProgram Edge_ErrorShader;
rts_glShaderProgram Node_ErrorShader;
rts_glShaderProgram Smooth_Shader;

//display lists
GLuint GT_FibersList=0;
GLuint T_FibersList=0;
GLuint GT_EndCaps=0;
GLuint T_EndCaps=0;
GLuint GT_NodesList=0;
GLuint T_NodesList=0;
GLuint T_PathList=0;
GLuint GT_PathList=0;

//drawing variables
int tube_subdivisions = 20;
float node_radius_factor = 0.7;
float fiber_radius_factor = 0.5;
float cull_test_case_threshold = 1.0;

#define DISPLAY_GT_NETWORK	1
#define DISPLAY_T_NETWORK	2
#define DISPLAY_GT_GRAPH	3
#define DISPLAY_T_GRAPH		4
#define DISPLAY_GT_SELECTED	5
#define DISPLAY_T_SELECTED	6

//menu options
#define NETCOMP_EXIT		0
#define DISPLAY_NETWORK		1
#define DISPLAY_GRAPH		2
#define DISPLAY_CONNECTED	3
#define DISPLAY_SELECTED	4
#define COLORMAP_ISOLUMINANT	5
#define COLORMAP_BLACKBODY		6
#define COLORMAP_BREWER			7
#define COLORMAP_POLAR_CIELAB	8
#define COLORMAP_RAINBOW		9
#define CULL_TEST_CASE			10
#define RECOMPUTE_METRIC		11


//fibers to render in Graph Mode
//list<int> T_DisplayEdges;
//list<int> GT_DisplayEdges;

int DisplayMode = DISPLAY_NETWORK;
float L0_pos[3];
float L1_pos[3];
//GLuint texColorMap=0;
rts_glTextureMap texColorMap;




void makeColormap(int ColorMapType = COLORMAP_BREWER)
{
	//if(texColorMap != 0)
	//	glDeleteTextures(1, &texColorMap);

	point3D<float>* ctrlPts;
	int num_points = 0;
	if(ColorMapType == COLORMAP_ISOLUMINANT)
	{
		//allocate memory for the colormap
		num_points = 2;
		ctrlPts = new point3D<float>[num_points];
		//memset(ctrlPts, 0, num_points*sizeof(point3D<float>));

		ctrlPts[0] = point3D<float>(0.0, 1.0, 0.0);
		ctrlPts[1] = point3D<float>(1.0, 0.0, 0.0);
	}
	else if(ColorMapType == COLORMAP_RAINBOW)
	{
		//allocate memory for the colormap
		num_points = 5;
		ctrlPts = new point3D<float>[num_points];
		//memset(ctrlPts, 0, num_points*sizeof(point3D<float>));

		//ctrlPts[0] = point3D<float>(0.7, 0, 0.7);
		ctrlPts[0] = point3D<float>(0, 0, 1);
		ctrlPts[1] = point3D<float>(0, 0.7, 0.7);
		ctrlPts[2] = point3D<float>(0, 1, 0);
		ctrlPts[3] = point3D<float>(0.7, 0.7, 0);
		ctrlPts[4] = point3D<float>(1, 0, 0);
	}
	else if(ColorMapType == COLORMAP_BLACKBODY)
	{
		//allocate memory for the colormap
		num_points = 4;
		ctrlPts = new point3D<float>[num_points];
		//memset(ctrlPts, 0, num_points*sizeof(point3D<float>));

		ctrlPts[0] = point3D<float>(0.0, 0.0, 0.0);
		ctrlPts[1] = point3D<float>(1.0, 0.0, 0.0);
		ctrlPts[2] = point3D<float>(1.0, 1.0, 0.0);
		ctrlPts[3] = point3D<float>(1.0, 1.0, 1.0);
	}
	else if(ColorMapType == COLORMAP_BREWER)
	{
		//allocate memory for the colormap
		num_points = 11;
		ctrlPts = new point3D<float>[num_points];
		//memset(ctrlPts, 0, num_points*sizeof(point3D<float>));

		ctrlPts[0] = point3D<float>(0.192157, 0.211765, 0.584314);
		ctrlPts[1] = point3D<float>(0.270588, 0.458824, 0.705882);
		ctrlPts[2] = point3D<float>(0.454902, 0.678431, 0.819608);
		ctrlPts[3] = point3D<float>(0.670588, 0.85098, 0.913725);
		ctrlPts[4] = point3D<float>(0.878431, 0.952941, 0.972549);
		ctrlPts[5] = point3D<float>(1, 1, 0.74902);
		ctrlPts[6] = point3D<float>(0.996078, 0.878431, 0.564706);
		ctrlPts[7] = point3D<float>(0.992157, 0.682353, 0.380392);
		ctrlPts[8] = point3D<float>(0.956863, 0.427451, 0.262745);
		ctrlPts[9] = point3D<float>(0.843137, 0.188235, 0.152941);
		ctrlPts[10] = point3D<float>(0.647059, 0, 0.14902);

	}
		else if(ColorMapType == COLORMAP_POLAR_CIELAB)
	{
		//allocate memory for the colormap
		num_points = 33;
		ctrlPts = new point3D<float>[num_points];
		//memset(ctrlPts, 0, num_points*sizeof(point3D<float>));

		ctrlPts[0] = point3D<float>(0.07514311, 0.468049805,1);
		ctrlPts[1] = point3D<float>(0.247872569, 0.498782363,1);
		ctrlPts[2] = point3D<float>(0.339526309, 0.528909511,1);
		ctrlPts[3] = point3D<float>(0.409505078, 0.558608486,1);
		ctrlPts[4] = point3D<float>(0.468487184, 0.588057293,1);
		ctrlPts[5] = point3D<float>(0.520796675, 0.617435078,1);
		ctrlPts[6] = point3D<float>(0.568724526, 0.646924167,1);
		ctrlPts[7] = point3D<float>(0.613686735, 0.676713218,1);
		ctrlPts[8] = point3D<float>(0.656658579, 0.707001303,1);
		ctrlPts[9] = point3D<float>(0.698372844, 0.738002964,1);
		ctrlPts[10] = point3D<float>(0.739424025, 0.769954435,1);
		ctrlPts[11] = point3D<float>(0.780330104, 0.803121429,1);
		ctrlPts[12] = point3D<float>(0.821573924, 0.837809045,1);
		ctrlPts[13] = point3D<float>(0.863634967, 0.874374691,1);
		ctrlPts[14] = point3D<float>(0.907017747, 0.913245283,1);
		ctrlPts[15] = point3D<float>(0.936129275, 0.938743558, 0.983038586);
		ctrlPts[16] = point3D<float>(0.943467973, 0.943498599, 0.943398095);
		ctrlPts[17] = point3D<float>(0.990146732, 0.928791426, 0.917447482);
		ctrlPts[18] = point3D<float>(1, 0.88332677, 0.861943246);
		ctrlPts[19] = point3D<float>(1, 0.833985467, 0.803839606);
		ctrlPts[20] = point3D<float>(1, 0.788626485, 0.750707739);
		ctrlPts[21] = point3D<float>(1, 0.746206642, 0.701389973);
		ctrlPts[22] = point3D<float>(1, 0.70590052, 0.654994046);
		ctrlPts[23] = point3D<float>(1, 0.667019783, 0.610806959);
		ctrlPts[24] = point3D<float>(1, 0.6289553, 0.568237474);
		ctrlPts[25] = point3D<float>(1, 0.591130233, 0.526775617);
		ctrlPts[26] = point3D<float>(1, 0.552955184, 0.485962266);
		ctrlPts[27] = point3D<float>(1, 0.513776083, 0.445364274);
		ctrlPts[28] = point3D<float>(1, 0.472800903, 0.404551679);
		ctrlPts[29] = point3D<float>(1, 0.428977855, 0.363073592);
		ctrlPts[30] = point3D<float>(1, 0.380759558, 0.320428137);
		ctrlPts[31] = point3D<float>(0.961891484, 0.313155629, 0.265499262);
		ctrlPts[32] = point3D<float>(0.916482116, 0.236630659, 0.209939162);

	}

	texColorMap.Init(ctrlPts, GL_TEXTURE_1D, num_points, 0, 0, GL_RGB, GL_RGB, GL_FLOAT);

	//glGenTextures(1, &texColorMap);
	//glBindTexture(GL_TEXTURE_1D, texColorMap);
	//glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, num_points, 0, GL_RGB, GL_FLOAT, ctrlPts);
}



point3D<float> HSLtoRGB(point3D<float> HSL)
{
	float H = HSL.x;
	float S = HSL.y;
	float L = HSL.z;
	
	float temp2;
	if(L < 0.5)
		temp2 = L*(1.0+S);
	else
		temp2 = L+S - L*S;
		
	float temp1 = 2.0*L - temp2;
	
	point3D<float> temp3(H+1.0/3.0, H, H-1.0/3.0);
	if(temp3.x < 0.0) temp3.x = temp3.x+1.0;
	if(temp3.y < 0.0) temp3.y = temp3.y+1.0;
	if(temp3.z < 0.0) temp3.z = temp3.z+1.0;
	
	if(temp3.x > 1.0) temp3.x = temp3.x - 1.0;
	if(temp3.y > 1.0) temp3.y = temp3.y - 1.0;
	if(temp3.z > 1.0) temp3.z = temp3.z - 1.0;
	
	point3D<float> result;
	if(6.0*temp3.x < 1.0) result.x = temp1 +(temp2 - temp1)*6.0*temp3.x;
	else if(2.0*temp3.x < 1.0) result.x = temp2;
	else if(3.0*temp3.x < 2.0) result.x = temp1+(temp2-temp1)*((2.0/3.0) - temp3.x)*6.0;
	else result.x = temp1;
	
	if(6.0*temp3.y < 1.0) result.y = temp1 +(temp2 - temp1)*6.0*temp3.y;
	else if(2.0*temp3.y < 1.0) result.y = temp2;
	else if(3.0*temp3.y < 2.0) result.y = temp1+(temp2-temp1)*((2.0/3.0) - temp3.y)*6.0;
	else result.y = temp1;
	
	if(6.0*temp3.z < 1.0) result.z = temp1 +(temp2 - temp1)*6.0*temp3.z;
	else if(2.0*temp3.z < 1.0) result.z = temp2;
	else if(3.0*temp3.z < 2.0) result.z = temp1+(temp2-temp1)*((2.0/3.0) - temp3.z)*6.0;
	else result.z = temp1;
	
	//result.a = 0.0;
	return result;
}
void ColorFibers()
{
	//srand(time(NULL));
	sequenceColors.clear();
	//for each fiber
	for(CoreGraphList::iterator i = coreGraph.begin(); i!=coreGraph.end(); i++)
	{
		float random_hue = (double)rand()/(double)RAND_MAX;
		//cout<<"Random Hue: "<<random_hue<<endl;
		float random_saturation = 1.0;//(double)rand()/(double)RAND_MAX;
		point3D<float> rgb = HSLtoRGB(point3D<float>(random_hue, random_saturation, 0.5));
		//point3D<float> rgb((double)rand()/(double)RAND_MAX, (double)rand()/(double)RAND_MAX, (double)rand()/(double)RAND_MAX);
		sequenceColors.push_back(rgb);
	}
	

}
void CenterCameraToSelected()
{
	if(coreGraph.size() == 0)
		return;

	//center the fiber in both networks
	point3D<float> min_pt(9999, 9999, 9999);
	point3D<float> max_pt(-9999, -9999, -9999);

	//iterate through the first edge sequence
	EdgeSequence::iterator i;
	int node;
	point3D<float> test;
	for(i=coreGraph[current_sequence].first.begin(); i!=coreGraph[current_sequence].first.end(); i++)
	{
		node = testNetwork->FiberList[*i].n0;
		test = testNetwork->NodeList[node].p;
		min_pt.x = min(test.x, min_pt.x);
		min_pt.y = min(test.y, min_pt.y);
		min_pt.z = min(test.z, min_pt.z);	
		max_pt.x = max(test.x, max_pt.x);
		max_pt.y = max(test.y, max_pt.y);
		max_pt.z = max(test.z, max_pt.z);

		node = testNetwork->FiberList[*i].n1;
		test = testNetwork->NodeList[node].p;
		min_pt.x = min(test.x, min_pt.x);
		min_pt.y = min(test.y, min_pt.y);
		min_pt.z = min(test.z, min_pt.z);	
		max_pt.x = max(test.x, max_pt.x);
		max_pt.y = max(test.y, max_pt.y);
		max_pt.z = max(test.z, max_pt.z);
	}
	point3D<float> middle = min_pt+0.5*(max_pt - min_pt);

	rts_glut_camera.LookAt(middle);

}
void IncrementSelectedFiber(int i)
{
	//get the currently selected fiber id
	if(coreGraph.size() <= 0)
		return;

	//get the number of fibers
	int end_id = coreGraph.size();

	current_sequence+=i;
	if(current_sequence >= end_id)
		current_sequence = 0;
	if(current_sequence < 0)
		current_sequence = coreGraph.size()-1;

	//print the selected edges
	EdgeSequence::iterator EdgeI;

	for(EdgeI = coreGraph[current_sequence].first.begin(); EdgeI != coreGraph[current_sequence].first.end(); EdgeI++)
		cout<<*EdgeI<<" ";
	cout<<"--->";
	for(EdgeI = coreGraph[current_sequence].second.begin(); EdgeI != coreGraph[current_sequence].second.end(); EdgeI++)
		cout<<*EdgeI<<" ";
	cout<<endl;


	CenterCameraToSelected();

}
void DrawNodeSphere(rtsFiberNetwork* network, int n, float radius)
{
	GLUquadricObj* quadric = gluNewQuadric();
	gluQuadricNormals(quadric, GLU_SMOOTH);

	glMatrixMode(GL_MODELVIEW);

	glPushMatrix();

	point3D<float> p;
	//glColor3f(network->NodeList[n].error, 0.0, 0.0);
	p = network->NodeList[n].p;

	glTranslatef(p.x, p.y, p.z);
	//glutSolidSphere(node_radius*standard_deviation, 20, 20);
	glTexCoord1f(network->NodeList[n].error);
	if(network->NodeList[n].color < 0)
		glColor3f(1.0, 0.0, 0.0);
	else
		glColor3f(1.0, 1.0, 1.0);
	gluSphere(quadric,radius,32,32);

	glPopMatrix();

}
void DrawNodeSpheres(rtsFiberNetwork* network, float radius)
{


	unsigned int n;
	for(n=0; n != network->FiberList.size(); n++)
	{
		if(!network->isCulled(n))
		{
			DrawNodeSphere(network, network->FiberList[n].n0,radius);
			DrawNodeSphere(network, network->FiberList[n].n1,radius);
		}
	}
}

void FrenetFrame(vector3D<float> &x, vector3D<float> &y, vector3D<float> &z)
{
	x = vector3D<float>(0.0, 0.0, 1.0);
	y = x.X(z);
	x = z.X(y);
	x.Normalize();
	y.Normalize();
	z.Normalize();
}

vector3D<float> GetColor(float error)
{
	//This function converts an error value to a color
	//The conversion is done by creating an HSV color from the error value and converting that HSV color to RGB
	float H = (240.0/60.0)*(1.0 - error);
	float S = 1.0;
	float V = 1.0;

	int i = floor(H);
	float f = H - i;
	if(i%2 == 0)
		f = 1-f;
	float m = V*(1 - S);
	float n = V*(1-S*f);
	switch(i)
	{
	case 0:
		return vector3D<float>(V, n, m);
	case 1:
		return vector3D<float>(n, V, m);
	case 2:
		return vector3D<float>(m, V, n);
	case 3:
		return vector3D<float>(m, n, V);
	case 4:
		return vector3D<float>(n, m, V);
	case 5:
		return vector3D<float>(V, m, n);
	default:
		return vector3D<float>(0, 0, 0);
	}

}

void DrawTube(point3D<float> p0, vector3D<float> d0, point3D<float> p1, vector3D<float> d1, float error0, float error1, float radius, int subdiv)
{
	
	//draw the first circle
	vector3D<float> x0, y0, z0, x1, y1, z1;
	z0 = d0;
	FrenetFrame(x0, y0, z0);

	z1 = d1;
	FrenetFrame(x1, y1, z1);

	float t_step = (2*3.14159)/subdiv;

	
	float u, v;

	//get the RGB color
	point3D<float> circle0, circle1;
	vector3D<float> RGB0, RGB1;
	vector3D<float> normal;
	
	//RGB0 = GetColor(color0);	
	//RGB1 = GetColor(color1);

	glBegin(GL_TRIANGLE_STRIP);
	for(int t=0; t<=subdiv; t++)
	{
		u = radius * cos(t*t_step);
		v = radius * sin(t*t_step);
		normal = u*x0 + v*y0;		
		circle0 = p0 + normal;
		normal.Normalize();

		glTexCoord1f(error0);
		//glColor4f(error0, 0.0, 0.0, 1.0);
		glNormal3f(normal.x, normal.y, normal.z);
		glVertex3f(circle0.x, circle0.y, circle0.z);

		normal = u*x1 + v*y1;
		circle1 = p1 + normal;
		normal.Normalize();

		glTexCoord1f(error1);
		//glColor4f(error1, 0.0, 0.0, 1.0);
		glNormal3f(normal.x, normal.y, normal.z);
		glVertex3f(circle1.x, circle1.y, circle1.z);

	}
	glEnd();
	CHECK_OPENGL_ERROR
}

void ExtrudeFiber(rtsFiberNetwork* network, int fiber, float radius)
{
	vector3D<float> x, y, z;
	point3D<float> p0, p1, p2, p3;
	vector3D<float> d1, d2;
	float e1, e2;

	//get the first point
	int node = network->FiberList[fiber].n0;
	p1 = network->NodeList[node].p;
	e1 = network->NodeList[node].error;

	//for each vertex in the fiber
	int num_points = (int)network->FiberList[fiber].pointList.size();
	for(int v=0; v<num_points; v++)
	{
		//get the next point
		p2 = network->FiberList[fiber].pointList[v];
		e2 = network->FiberList[fiber].errorList[v];
		
		if(v<num_points-1)
			p3 = network->FiberList[fiber].pointList[v+1];
		else
		{
			node = network->FiberList[fiber].n1;
			p3 = network->NodeList[node].p;
		}

		d2 = p3-p1;

		//compute the fiber derivatives at p1 and p2
		if(v==0)	//if this is the first fiber
			d1 = p2 - p1;
		else
		{
			d1 = p2 - p0;
		}

		DrawTube(p1, d1, p2, d2, e1, e2, radius, tube_subdivisions);

		//shift
		p0 = p1;
		p1 = p2;
		e1 = e2;
	}
	//make the last tube

	//if there were any points in the pointlist
	if(num_points > 0)
	{
		p2 = p3;
		node = network->FiberList[fiber].n1;
		e2 = network->NodeList[node].error;
		d1 = p2-p0;
		d2 = p2-p1;
		DrawTube(p1, d1, p2, d2, e1, e2, radius, tube_subdivisions);
	}
	//if there are only the two node points
	else
	{
		node = network->FiberList[fiber].n1;
		p2 = network->NodeList[node].p;
		e2 = network->NodeList[node].error;
		d1 = p2 - p1;
		d2 = p2 - p1;
		DrawTube(p1, d1, p2, d2, e1, e2, radius, tube_subdivisions);
	}
}

void DrawLineFiber(rtsFiberNetwork* network, int f)
{
	point3D<float> p;
	int node = network->FiberList[f].n0;
	p = network->NodeList[node].p;

	glBegin(GL_LINE_STRIP);
	glVertex3f(p.x, p.y, p.z);
	for(int v=0; v!=network->FiberList[f].pointList.size(); v++)
	{
		p = network->FiberList[f].pointList[v];
		glVertex3f(p.x, p.y, p.z);		
	}
	node = network->FiberList[f].n1;
	p = network->NodeList[node].p;
	glVertex3f(p.x, p.y, p.z);
	glEnd();

}
void DrawLineNetwork(rtsFiberNetwork* network)
{

	int num_fibers = network->FiberList.size();
	for(int f = 0; f < num_fibers; f++)
	{
		/*if(network->FiberList[f].mapped_to == -1)
			glColor3f(0.0,0.0, 0.0);
		else
			glColor3f(1.0, 0.0, 0.0);
		*/
		DrawLineFiber(network, f);
		CHECK_OPENGL_ERROR
	}


}
void DrawGraphNodes(rtsFiberNetwork* network)
{
	//renders graph nodes, colored based on their node color
	glMatrixMode(GL_MODELVIEW);

	unsigned int n;
	for(n=0; n != network->NodeList.size(); n++)
	{
		glPushMatrix();

		point3D<float> p;
		/*if(network->NodeList[n].color < 0)
			glColor3f(1.0, 0.0, 0.0);
		else
			glColor3f(0.0, 1.0, 0.0);
			*/
		//glColor3f(network->NodeList[n].error, 0.0, 0.0);
		p = network->NodeList[n].p;

		glTranslatef(p.x, p.y, p.z);
		glutSolidSphere(node_radius_factor*sigmaC, 20, 20);

		glPopMatrix();
	}
}
void DrawFiberSequence(rtsFiberNetwork* network, EdgeSequence sequence, float fiber_radius, float node_radius)
{
	//glClear(GL_DEPTH_BUFFER_BIT);
	for(EdgeSequence::iterator i = sequence.begin(); i != sequence.end(); i++)
	{
		ExtrudeFiber(network, *i, fiber_radius);
		glPushAttrib(GL_CURRENT_BIT);
		DrawNodeSphere(network, network->FiberList[*i].n0, node_radius);
		DrawNodeSphere(network, network->FiberList[*i].n1, node_radius);
		glPopAttrib();
	}


}
GLuint CreateFiberDisplayList(rtsFiberNetwork* network, float radius)
{
	GLuint result = glGenLists(1);
	glNewList(result, GL_COMPILE);

	int num_fibers = network->FiberList.size();
	for(int f = 0; f < num_fibers; f++)
	{
		if(!network->isCulled(f))
		{
			ExtrudeFiber(network, f, radius);
			CHECK_OPENGL_ERROR
		}
	}
	glEndList();
	return result;
}

GLuint CreateNodeDisplayList(rtsFiberNetwork* network, float radius)
{
	GLuint result = glGenLists(1);
	glNewList(result, GL_COMPILE);
	DrawNodeSpheres(network, radius);
	glEndList();
	return result;
}


void CreateFiberPathLists(float fiber_radius, float node_radius)
{
	GT_PathList = glGenLists(1);
	glNewList(GT_PathList, GL_COMPILE);

	if(coreGraph.size() > 0)
	{
		for(unsigned int i=0; i<sequenceColors.size(); i++)
		{
			point3D<float> rgb = sequenceColors[i];
			glColor3f(rgb.x, rgb.y, rgb.z);
			DrawFiberSequence(goldNetwork, coreGraph[i].second, fiber_radius, node_radius);
		}
	}
	glEndList();

	T_PathList = glGenLists(1);
	glNewList(T_PathList, GL_COMPILE);

	if(coreGraph.size() > 0)
	{
		for(unsigned int i=0; i<sequenceColors.size(); i++)
		{
			point3D<float> rgb = sequenceColors[i];
			glColor3f(rgb.x, rgb.y, rgb.z);
			DrawFiberSequence(testNetwork, coreGraph[i].first, fiber_radius, node_radius);
		}
	}
	glEndList();

}
void CreateDisplayLists()
{
	if(GT_FibersList != 0)
		glDeleteLists(GT_FibersList, 1);
	if(T_FibersList != 0)
		glDeleteLists(T_FibersList, 1);
	if(GT_NodesList != 0)
		glDeleteLists(GT_NodesList, 1);
	if(T_NodesList != 0)
		glDeleteLists(T_NodesList, 1);
	if(GT_EndCaps != 0)
		glDeleteLists(GT_EndCaps, 1);
	if(T_EndCaps != 0)
		glDeleteLists(T_EndCaps, 1);

	//create the display lists
	GT_FibersList = CreateFiberDisplayList(goldNetwork, sigmaG*fiber_radius_factor);
	T_FibersList = CreateFiberDisplayList(testNetwork, sigmaG*fiber_radius_factor);
	GT_NodesList = CreateNodeDisplayList(goldNetwork, sigmaG*node_radius_factor);
	T_NodesList = CreateNodeDisplayList(testNetwork, sigmaG*node_radius_factor);
	GT_EndCaps = CreateNodeDisplayList(goldNetwork, sigmaG*fiber_radius_factor);
	T_EndCaps = CreateNodeDisplayList(testNetwork, sigmaG*fiber_radius_factor);

	if(GT_PathList != 0)
		glDeleteLists(GT_PathList,1);
	if(T_PathList != 0)
		glDeleteLists(T_PathList,1);
	CreateFiberPathLists(fiber_radius_factor*sigmaG, node_radius_factor*sigmaG);
}
void RenderViewCamera()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//compute the aspect ratio
	float aspect_ratio = (float)glutGet(GLUT_WINDOW_WIDTH)/2.0/(float)glutGet(GLUT_WINDOW_HEIGHT);
	gluPerspective(rts_glut_camera.getFOV(), aspect_ratio, network_span/10, network_span*10);

	//render the camera
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	point3D<float> camera_position = rts_glut_camera.getPosition();
	vector3D<float> camera_up = rts_glut_camera.getUp();
	point3D<float> camera_lookat = rts_glut_camera.getLookAt();
	gluLookAt(camera_position.x,
			  camera_position.y,
			  camera_position.z,
			  camera_lookat.x,
			  camera_lookat.y,
			  camera_lookat.z,
			  camera_up.x,
			  camera_up.y,
			  camera_up.z);


	//get the light positions (lights move with the camera)
	vector3D<float> up = rts_glut_camera.getUp();
	vector3D<float> dir = rts_glut_camera.getDirection();
	vector3D<float> side = dir.X(up);
	L0_pos[0] = side.x;
	L0_pos[1] = side.y;
	L0_pos[2] = side.z;

	L1_pos[0] = -dir.x;
	L1_pos[1] = -dir.y;
	L1_pos[2] = -dir.z;

	//scale the viewport to the network
	/*vector3D<float> span = max_point - min_point;
	float scale = 1.0/span.Length();
	//compute center point
	point3D<float> center = min_point + 0.5*span;
	glScalef(scale, scale, scale);
	glTranslatef(-center.x, -center.y, -center.z);*/


}
void RecenterCamera()
{
	point3D<float> min_point0 = goldNetwork->min_pos;
	point3D<float> min_point1 = testNetwork->min_pos;

	point3D<float> max_point0 = goldNetwork->max_pos;
	point3D<float> max_point1 = testNetwork->max_pos;

	point3D<float> min_point(min(min_point0.x, min_point1.x), min(min_point0.y, min_point1.y), min(min_point0.z, min_point1.z));
	point3D<float> max_point(max(max_point0.x, max_point1.x), max(max_point0.y, max_point1.y), max(max_point0.z, max_point1.z));
	point3D<float> center_point = min_point + 0.5*(max_point - min_point);


	network_span = (max_point - point3D<float>(0, 0, 0)).Length();
	rts_glut_camera.setPosition(0, 0, 3*network_span);
	rts_glut_camera.LookAt(center_point, vector3D<float>(0.0, 1.0, 0.0));

}
void DrawNetwork(int Display)
{
	//draw the network fibers
	switch(Display)
	{
	case DISPLAY_GT_NETWORK:
		Edge_ErrorShader.UpdateGlobalUniforms();
		Edge_ErrorShader.BeginProgram();
		glCallList(GT_FibersList);
		Edge_ErrorShader.EndProgram();

		Node_ErrorShader.UpdateGlobalUniforms();
		Node_ErrorShader.BeginProgram();
		glCallList(GT_EndCaps);
		Node_ErrorShader.EndProgram();
		break;
	case DISPLAY_T_NETWORK:
		Edge_ErrorShader.UpdateGlobalUniforms();
		Edge_ErrorShader.BeginProgram();
		glCallList(T_FibersList);	
		Edge_ErrorShader.EndProgram();

		Node_ErrorShader.UpdateGlobalUniforms();
		Node_ErrorShader.BeginProgram();
		glCallList(T_EndCaps);	
		Node_ErrorShader.EndProgram();
		break;
	case DISPLAY_GT_GRAPH:
		glColor3f(1.0, 1.0, 1.0);
		Smooth_Shader.UpdateGlobalUniforms();
		Smooth_Shader.BeginProgram();
		glCallList(GT_PathList);
		glColor3f(1.0, 1.0, 1.0);
		glCallList(GT_FibersList);
		glCallList(GT_NodesList);
		Smooth_Shader.EndProgram();
		break;
	case DISPLAY_T_GRAPH:
		glColor3f(1.0, 1.0, 1.0);
		Smooth_Shader.UpdateGlobalUniforms();
		Smooth_Shader.BeginProgram();
		glCallList(T_PathList);
		glColor3f(1.0, 1.0, 1.0);
		glCallList(T_FibersList);
		glCallList(T_NodesList);
		Smooth_Shader.EndProgram();
		break;
	case DISPLAY_GT_SELECTED:
		glColor3f(1.0, 1.0, 1.0);
		Smooth_Shader.UpdateGlobalUniforms();
		Smooth_Shader.BeginProgram();
		glCallList(GT_FibersList);
		glCallList(GT_NodesList);

		glClear(GL_DEPTH_BUFFER_BIT);
		glColor3f(1.0, 0.0, 1.0);
		DrawFiberSequence(goldNetwork, coreGraph[current_sequence].second, fiber_radius_factor*sigmaG, node_radius_factor*sigmaG);
		Smooth_Shader.EndProgram();
		break;
	case DISPLAY_T_SELECTED:
		glColor3f(1.0, 1.0, 1.0);
		Smooth_Shader.UpdateGlobalUniforms();
		Smooth_Shader.BeginProgram();
		glCallList(T_FibersList);
		glCallList(T_NodesList);

		glClear(GL_DEPTH_BUFFER_BIT);
		glColor3f(1.0, 0.0, 1.0);
		DrawFiberSequence(testNetwork, coreGraph[current_sequence].first, fiber_radius_factor*sigmaG, node_radius_factor*sigmaG);
		Smooth_Shader.EndProgram();
		break;
	default:
		break;
	}

}




void MyDisplayFunction()
{
	//glClearColor(0.95, 1.0, 0.95, 1.0);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);

	//set left viewport
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH)/2, glutGet(GLUT_WINDOW_HEIGHT));

	RenderViewCamera();

	
	if(DisplayMode == DISPLAY_NETWORK)
		DrawNetwork(DISPLAY_GT_NETWORK);
	else if(DisplayMode == DISPLAY_GRAPH)
		DrawNetwork(DISPLAY_GT_GRAPH);
	else if(DisplayMode == DISPLAY_SELECTED)
		DrawNetwork(DISPLAY_GT_SELECTED);

	//set right viewport
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(glutGet(GLUT_WINDOW_WIDTH)/2, 0, glutGet(GLUT_WINDOW_WIDTH)/2, glutGet(GLUT_WINDOW_HEIGHT));
	
	RenderViewCamera();
	//glClear(GL_COLOR_BUFFER_BIT);
	//glClear(GL_DEPTH_BUFFER_BIT);

	if(DisplayMode == DISPLAY_NETWORK)
		DrawNetwork(DISPLAY_T_NETWORK);
	else if(DisplayMode == DISPLAY_GRAPH)
		DrawNetwork(DISPLAY_T_GRAPH);
	else if(DisplayMode == DISPLAY_SELECTED)
		DrawNetwork(DISPLAY_T_SELECTED);


	
	//glutSwapBuffers();
}

void GlutMenuCallback(int option)
{
	if(option == NETCOMP_EXIT)
		exit(1);
	if(option >= DISPLAY_NETWORK && option <= DISPLAY_SELECTED)
		DisplayMode = option;
	if(option >= COLORMAP_ISOLUMINANT && option <= COLORMAP_RAINBOW)
		makeColormap(option);
	if(option == CULL_TEST_CASE)
	{
		//get the new threshold
		cout<<"Enter new test case fiber threshold [0 1]: ";
		float cull_value;
		cin>>cull_value;
		testNetwork->setCullValue(cull_value);

		//re-create the display lists
		ComputeNetMets();
		CreateDisplayLists();

	}
	if(option == RECOMPUTE_METRIC)
	{
		cout<<"Please enter a sigma value: ";
		cin>>sigmaG;
		sigmaC = sigmaG;

		ComputeNetMets();
		CreateDisplayLists();
	}

	

	/*
	switch(option)
	{
	case NETCOMP_EXIT:
		exit(1);
		break;
	case NETCOMP_VIEW_GT:
		DisplayMode = DISPLAY_GT_NETWORK;
		break;
	case NETCOMP_VIEW_T:
		DisplayMode = DISPLAY_T_NETWORK;
		break;
	default:
		break;
	}*/

}