#include <rtsPoint3d.h>
#include <vector>
#include <string>
//#include <Windows.h>

//#include "ImplicitCalculations.h"
#include "objJedi.h"
#include "DrawingFunctions.h"
#include "rtsGraph.h"
#include "rtsFilename.h"

#define PRINT_LOG	1


using namespace std;

rtsFiberNetwork* goldNetwork;
GLint goldNodeDL;

rtsFiberNetwork* testNetwork;

//standard deviation of the gaussian envelope surrounding the networks
float sigmaG, sigmaC;

int glutwindow_id;

void GlutMenuCallback(int option);

void LoadNetworks(string gold_filename, string test_filename)
{
	goldNetwork = new rtsFiberNetwork();
	testNetwork = new rtsFiberNetwork();

	cout<<endl;
	cout<<"Gold Standard:"<<endl;
	goldNetwork->LoadSWC(gold_filename);
	cout<<"Nodes: "<<goldNetwork->NodeList.size()<<endl;
	cout<<"Fibers: "<<goldNetwork->FiberList.size()<<endl;
	cout<<"Total Length: "<<goldNetwork->getTotalLength()<<endl;

	cout<<endl;

	cout<<"Test Network:"<<endl;
	testNetwork->LoadSWC(test_filename);
	cout<<"Nodes: "<<testNetwork->NodeList.size()<<endl;
	cout<<"Fibers: "<<testNetwork->FiberList.size()<<endl;
	cout<<"Total Length: "<<testNetwork->getTotalLength()<<endl;

	cout<<endl;

	//focusNetwork = testNetwork;
}

void SpecialKeys(int key, int x, int y)
{
	//rtsQuaternion<float> new_rotation;

	//if(key == GLUT_KEY_UP)
		//rts_glut_camera.OrbitFocus(0, d_angle);
	//if(key == GLUT_KEY_DOWN)
		//rts_glut_camera.OrbitFocus(0, -d_angle);
	if(key == GLUT_KEY_LEFT)
		IncrementSelectedFiber(-1);
		//rts_glut_camera.OrbitFocus(d_angle, 0.0);
	if(key == GLUT_KEY_RIGHT)
		IncrementSelectedFiber(1);
		//rts_glut_camera.OrbitFocus(-d_angle, 0.0);
}

void PrintHelp()
{
	cout<<"'.' = Center Camera"<<endl;
	cout<<"'z' = Reduce Radius"<<endl;
	cout<<"'x' = Increase Radius"<<endl;
	cout<<"'c' = Color Segments"<<endl;
}

void KeyboardFunction(unsigned char key, int x, int y)
{
	if(key == '.')
		RecenterCamera();
	if(key == '?')
		PrintHelp();
	if(key == 'z')
	{
		fiber_radius_factor += 0.1;
		node_radius_factor += 0.1;
		CreateDisplayLists();
	}
	if(key == 'x')
	{
		fiber_radius_factor -= 0.1;
		node_radius_factor -= 0.1;
		CreateDisplayLists();
	}
	if(key == 'c')
	{
		ColorFibers();
		CreateDisplayLists();
	}
	if(key == 'q')
	{
		exit(1);
	}

}

void ComputeNetMets()
{

	testNetwork->Resample(sigmaG);
	testNetwork->SubdivideNetwork(sigmaG);

	goldNetwork->Resample(sigmaG);
	goldNetwork->SubdivideNetwork(sigmaG);


	//compare the two networks and get the core graph
	coreGraph = goldNetwork->CompareNetworks(testNetwork, sigmaG, sigmaC);

	//color each of the core graph fiber sequences
	ColorFibers();

}

//#include <direct.h>
   #include <stdio.h>
   #include <errno.h>
int main(int argc, char* argv[])
{
	/*char pBuf[2048];
	GetModuleFileNameA(NULL, pBuf, 2048);
	string buffer = pBuf;
	rtsFilename filename = buffer;
	string ExecutableDirectory = filename.getDirectory();
	ExecutableDirectory+="\\";*/




	string gold_filename;
	string test_filename;
	sigmaG = 25;
	sigmaC = 25;
	int resolution = 512;
	float epsilon = 0.1;
	if(argc < 3)
	{
		gold_filename = "00_GT.obj"; test_filename = "00_T.obj"; sigmaG = sigmaC = 25.0;
		sigmaG = 25; sigmaC = 25.0;

	}
	else if(argc == 3)
	{
		gold_filename = argv[1];
		cout<<"Gold: "<<gold_filename<<endl;
		test_filename = argv[2];
		cout<<"Test: "<<test_filename<<endl;
		cout<<"Please enter a sigma value: ";
		cin>>sigmaG;
		sigmaC = sigmaG;
	}
	else if(argc == 4)
	{
		gold_filename = argv[1];
		cout<<"Gold: "<<gold_filename<<endl;
		test_filename = argv[2];
		cout<<"Test: "<<test_filename<<endl;
		sigmaG = sigmaC = atof(argv[3]);
	}

	//load the network files
	testNetwork = new rtsFiberNetwork();
	testNetwork->LoadFile(test_filename);

	goldNetwork = new rtsFiberNetwork();
	goldNetwork->LoadFile(gold_filename);




	//GLUT stuff
	//menus
	rts_glutInitialize("Network Comparison", 1000, 500);

	int mnuMain = glutCreateMenu(GlutMenuCallback);
	int mnuColormap = glutCreateMenu(GlutMenuCallback);

	glutSetMenu(mnuMain);
	glutAddMenuEntry("Network", DISPLAY_NETWORK);
	glutAddMenuEntry("Graph", DISPLAY_GRAPH);
	glutAddMenuEntry("Connected", DISPLAY_CONNECTED);
	glutAddMenuEntry("Selected", DISPLAY_SELECTED);

	//create the colormap menu
	glutAddSubMenu("Colormap", mnuColormap);
	glutSetMenu(mnuColormap);
	glutAddMenuEntry("Isoluminant", COLORMAP_ISOLUMINANT);
	glutAddMenuEntry("Black Body", COLORMAP_BLACKBODY);
	glutAddMenuEntry("Brewer", COLORMAP_BREWER);
	glutAddMenuEntry("Rainbow", COLORMAP_RAINBOW);

	glutSetMenu(mnuMain);
	glutAddMenuEntry("Exit", NETCOMP_EXIT);
	glutAddMenuEntry("Cull Test-Case", CULL_TEST_CASE);
	glutAddMenuEntry("Recompute NetMets", RECOMPUTE_METRIC);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	//keyboard stuff
	glutSpecialFunc(SpecialKeys);
	glutKeyboardFunc(KeyboardFunction);

	//camera
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

	//load the shaders
	makeColormap();
	//create strings for all of the shader filenames
	string SmoothShaderVertexFilename = "";
		SmoothShaderVertexFilename += "SmoothShader_Vertex.glsl";
	string SmoothShaderFragmentFilename = "";
		SmoothShaderFragmentFilename += "SmoothShader_Fragment.glsl";
	string ErrorMapFragmentFilename = "";
		ErrorMapFragmentFilename += "ErrorMap_Fragment.glsl";

	//Edge_ErrorShader.Init();
	//ErrorShader.Clean();
	Edge_ErrorShader.AttachShader(GL_VERTEX_SHADER, SmoothShaderVertexFilename.c_str());
	//Edge_ErrorShader.AttachShader(GL_FRAGMENT_SHADER, "ErrorShader_Fragment.glsl");
	Edge_ErrorShader.AttachShader(GL_FRAGMENT_SHADER, ErrorMapFragmentFilename.c_str());
	Edge_ErrorShader.Compile();
	//cout<<"Edge Error Shader Log-------------------------"<<endl;
	Edge_ErrorShader.PrintLog();
	Edge_ErrorShader.Link();
	Edge_ErrorShader.PrintLog();
	//Edge_ErrorShader.PrintUniforms();
	Edge_ErrorShader.AttachGlobalUniform("L0_pos", L0_pos);
	Edge_ErrorShader.AttachGlobalUniform("L1_pos", L1_pos);
	Edge_ErrorShader.AttachTextureMap("colorMap", texColorMap);
	//Edge_ErrorShader.UpdateGlobalUniforms();

	Node_ErrorShader.AttachShader(GL_VERTEX_SHADER, SmoothShaderVertexFilename.c_str());
	Node_ErrorShader.AttachShader(GL_FRAGMENT_SHADER, ErrorMapFragmentFilename.c_str());
	Node_ErrorShader.Compile();
	//cout<<"Node Error Shader Log-------------------------"<<endl;
	Node_ErrorShader.PrintLog();
	Node_ErrorShader.Link();
	Node_ErrorShader.PrintLog();
	Node_ErrorShader.AttachGlobalUniform("L0_pos", L0_pos);
	Node_ErrorShader.AttachGlobalUniform("L1_pos", L1_pos);
	Node_ErrorShader.AttachTextureMap("colorMap", texColorMap);

	Smooth_Shader.AttachShader(GL_VERTEX_SHADER, SmoothShaderVertexFilename.c_str());
	Smooth_Shader.AttachShader(GL_FRAGMENT_SHADER, SmoothShaderFragmentFilename.c_str());
	Smooth_Shader.Compile();
	//cout<<"Color Shader Log-------------------------"<<endl;
	Smooth_Shader.PrintLog();
	Smooth_Shader.Link();
	Smooth_Shader.PrintLog();
	Smooth_Shader.AttachGlobalUniform("L0_pos", L0_pos);
	Smooth_Shader.AttachGlobalUniform("L1_pos", L1_pos);

	//compute the metric
	ComputeNetMets();

	CreateDisplayLists();



	rts_glutStart(MyDisplayFunction);


	cout<<"Press <Enter> to Exit..."<<endl;
	cin.get();
	cin.get();

	return 1;

}
