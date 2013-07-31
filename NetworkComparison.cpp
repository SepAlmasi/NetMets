#include <rtsPoint3d.h>
#include <vector>
#include <string>

//AnyOption library for parsing command-line options:
//http://www.hackorama.com/anyoption/
#include "anyoption.h"

//#include "ImplicitCalculations.h"
#include "objJedi.h"
#include "DrawingFunctions.h"
#include "rtsGraph.h"
#include "rtsFilename.h"

#define PRINT_LOG	1


using namespace std;

bool displayGui = false;
bool displayVerbose = false;

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
	if(key == GLUT_KEY_LEFT)
		IncrementSelectedFiber(-1);
	if(key == GLUT_KEY_RIGHT)
		IncrementSelectedFiber(1);
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

void InitializeOpenGL()
{
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
	Edge_ErrorShader.AttachGlobalUniform("L1_pos", L1_pos);
	Edge_ErrorShader.AttachTextureMap("colorMap", texColorMap);


	Node_ErrorShader.AttachShader(GL_VERTEX_SHADER, SmoothShaderVertexFilename.c_str());
	Node_ErrorShader.AttachShader(GL_FRAGMENT_SHADER, ErrorMapFragmentFilename.c_str());
	Node_ErrorShader.Compile();
	Node_ErrorShader.PrintLog();
	Node_ErrorShader.Link();
	Node_ErrorShader.PrintLog();
	Node_ErrorShader.AttachGlobalUniform("L1_pos", L1_pos);
	Node_ErrorShader.AttachTextureMap("colorMap", texColorMap);

	Smooth_Shader.AttachShader(GL_VERTEX_SHADER, SmoothShaderVertexFilename.c_str());
	Smooth_Shader.AttachShader(GL_FRAGMENT_SHADER, SmoothShaderFragmentFilename.c_str());
	Smooth_Shader.Compile();
	Smooth_Shader.PrintLog();
	Smooth_Shader.Link();
	Smooth_Shader.PrintLog();
	Smooth_Shader.AttachGlobalUniform("L0_pos", L0_pos);
	Smooth_Shader.AttachGlobalUniform("L1_pos", L1_pos);

}

void ComputeNetMets()
{

	testNetwork->Resample(sigmaG);
	testNetwork->SubdivideNetwork(sigmaG);

	goldNetwork->Resample(sigmaG);
	goldNetwork->SubdivideNetwork(sigmaG);


	//compare the two networks and get the core graph
	float gFPR, gFNR, cFPR, cFNR;
	coreGraph = goldNetwork->CompareNetworks(testNetwork, sigmaG, sigmaC, gFPR, gFNR, cFPR, cFNR);

	//output error metrics
	if(displayVerbose)
	{
        cout<<"GEOMETRY++++++++++++"<<endl;
        cout<<"False Positive Rate: "<<gFPR<<endl;
        cout<<"False Negative Rate: "<<gFNR<<endl;
        cout<<"CONNECTIVITY++++++++"<<endl;
        cout<<"False Positive Rate: "<<cFPR<<endl;
        cout<<"False Negative Rate: "<<cFNR<<endl;
    }
    else
    {
        cout<<gFPR<<'\t'<<gFNR<<'\t'<<cFPR<<'\t'<<cFNR<<endl;
    }

	//color each of the core graph fiber sequences
	ColorFibers();

}

#include <stdio.h>
#include <errno.h>
int main(int argc, char* argv[])
{
    //Create an AnyOption object
    AnyOption* opt = new AnyOption();
    opt->addUsage( "" );
    opt->addUsage( "Usage: NetMets [OPTION]... [GROUND TRUTH FILE] [TEST CASE FILE]" );
    opt->addUsage( "" );
    opt->addUsage( " -h  --help         Prints this help " );
    opt->addUsage( " -s  --sigma        Input sigma value (default = 25)");
    opt->addUsage( " -g  --gui          Display NetMets GUI " );
    opt->addUsage( " -v  --verbose      Verbose output " );
    opt->addUsage( "" );

    opt->setFlag(  "help", 'h' );   /* a flag (takes no argument), supporting long and short form */
    opt->setFlag(  "gui", 'g' );
    opt->setOption("sigma", 's');
    opt->setFlag( "verbose", 'v');

    /* go through the command line and get the options  */
    opt->processCommandArgs( argc, argv );

    //display the help
    if( opt->getFlag( "help" ) || opt->getFlag( 'h' ) )
    {
                opt->printUsage();
                return 0;
    }

    //set a flag to display the GUI if requested
    if( opt->getFlag( "gui" ) || opt->getFlag( 'g' ) )
                displayGui = true;
    if( opt->getFlag("verbose") || opt->getFlag('v'))
                displayVerbose = true;

    //set the sigma value based on user input
    sigmaG = 25;
	sigmaC = 25;

	if( opt->getValue( 's' ) != NULL  || opt->getValue( "sigma" ) != NULL  )
	{
        sigmaG = atof(opt->getValue("sigma"));
        sigmaC = sigmaG;

    }
        	//cout << "size = " << atof(opt->getValue("sigma")) << endl ;

    //get the filename arguments
    int nArgs = opt->getArgc();
    char** sArgs = (char**)malloc(nArgs);
    for(int a=0; a<nArgs; a++)
        sArgs[a] = opt->getArgv(a);


	string gold_filename;
	string test_filename;

	int resolution = 512;
	float epsilon = 0.1;
	if(nArgs < 2)
	{
		gold_filename = "00_GT.obj"; test_filename = "00_T.obj"; sigmaG = sigmaC = 25.0;
		//sigmaG = 25; sigmaC = 25.0;

	}
	else
	{
        gold_filename = sArgs[0];
		test_filename = sArgs[1];
	}


	//load the network files
	testNetwork = new rtsFiberNetwork();
	testNetwork->LoadFile(test_filename);

	goldNetwork = new rtsFiberNetwork();
	goldNetwork->LoadFile(gold_filename);






	//compute the metric
	ComputeNetMets();



	//if the user does not want the GUI displayed, just return
    if(!displayGui) return 0;

    //initialize OpenGL
	InitializeOpenGL();

    //create display lists
	CreateDisplayLists();


    //set up GLUT and start
	rts_glutStart(MyDisplayFunction);

	return 0;

}
