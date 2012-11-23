#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <rtsLinearAlgebra.h>
//#include "DrawingFunctions.h"
#include "rtsFiberNetwork.h"
#include <octree\octree.h>
#include <rts_itkVolume.h>
#include <rtsDTGrid3D.h>
#include <rtsDTGridFastSweeping.h>
#include <algorithm>

/***************************************************
New Stuff
***************************************************/
//networks for analysis
rtsFiberNetwork* goldNetwork;
rtsFiberNetwork* testNetwork;
rtsFiberNetwork* focusNetwork;


//Implicit network representations
rtsDTGrid3D<float>* goldGrid;
rtsDTGrid3D<float>* testGrid;
rtsDTGrid3D<float>* focusGrid;

//lexicographic lists, contains LexIDs for each sampled point as well as identifiers to the associated fiber
struct LexID
{
	unsigned long ID;
	unsigned int fiber;

	bool operator<(LexID &rhs)
	{
		if(ID < rhs.ID)
			return true;
		return false;
	}
	bool operator==(LexID &rhs)
	{
		if(ID == rhs.ID)
			return true;
		return false;
	}
};
list<LexID> goldLexList;
list<LexID> testLexList;

//standard deviation of the implicit envelope
float STD;

//Size of each grid voxel.  This value is defined in Network Space.
float VOXEL_SIZE;

//Theoretical size of the grid.  Used to compute the lexicographic identifiers.
vector3D<int> GRID_DIM;

//Transformation matrix that converts Network Space coordinates into Grid Space coordinates.
matrix4x4<float> toGRID;

//Dilation parameter describes the radius of the envelope around the network (in voxels)
int DILATION;

void ComputeGridProperties()
{
	/*Computes the following grid properties:
	VOXEL_SIZE = the size (along one side) of each voxel "bin"
	GRID_DIM = the theoretical dimensions of the grid used to represent the network
	toGRID = transformation used to convert a point on the network into "grid space"
	*/

	//pick a voxel size based on the standard daviation of the envelope that captures the shape of the Gaussian
	//lets say that a voxel is 1 standard-deviation wide
	VOXEL_SIZE = STD*1.0;
	//dilate the DT Grid to 3 standard deviations
	DILATION = (3.0*STD)/VOXEL_SIZE;
	

	//find the minimum and maximum points in both networks to get the theoretical grid size
	//find the bounding box for both networks
	point3D<float> net_min;
	net_min.x = min(goldNetwork->min_pos.x, testNetwork->min_pos.x);
	net_min.y = min(goldNetwork->min_pos.y, testNetwork->min_pos.y);
	net_min.z = min(goldNetwork->min_pos.z, testNetwork->min_pos.z);

	point3D<float> net_max;
	net_max.x = max(goldNetwork->max_pos.x, testNetwork->max_pos.x);
	net_max.y = max(goldNetwork->max_pos.y, testNetwork->max_pos.y);
	net_max.z = max(goldNetwork->max_pos.z, testNetwork->max_pos.z);

	//compute the spatial size of the network
	vector3D<float> net_size = net_max - net_min;

	//compute the theoretical size of the grid
	GRID_DIM.x = ceil(net_size.x / VOXEL_SIZE);
	GRID_DIM.y = ceil(net_size.y / VOXEL_SIZE);
	GRID_DIM.z = ceil(net_size.z / VOXEL_SIZE);

	//compute the transformation from Network Space to Grid Space
	matrix4x4<float> scale;
	scale. SetIdentity();
	scale.SetScale(1.0/VOXEL_SIZE, 1.0/VOXEL_SIZE, 1.0/VOXEL_SIZE);

	matrix4x4<float> trans;
	trans.SetIdentity();
	trans.SetTranslation(-net_min.x, -net_min.y, -net_min.z);

	toGRID.SetIdentity();
	toGRID = scale*trans;

	cout<<"Voxel Size: "<<VOXEL_SIZE<<endl;
	cout<<"Grid Dimensions: "<<GRID_DIM.x<<","<<GRID_DIM.y<<","<<GRID_DIM.z<<endl;
}

void RasterizeSegmentLex(int f, point3D<float> p0, point3D<float> p1, list<LexID>* destList)
{
	/*Resamples a segment and stores the resulting points as
	lexicographic values in LexicographicList.
	*/

	//transform the Network Space coordinates to Grid Space
	point3D<float> grid_p0 = toGRID*p0;
	point3D<float> grid_p1 = toGRID*p1;

	//find the direction of travel
	vector3D<float> v = grid_p1 - grid_p0;

	//set the step size to the voxel size
	int length = (int)v.Length();
	v.Normalize();

	//start at p0, continue until p1
	point3D<float> p = grid_p0;

	int l;
	LexID lex;

	
	for(l=0; l<=length; l++)
	{
		lex.ID = (unsigned long)p.x*GRID_DIM.y*GRID_DIM.z + (unsigned long)p.y*GRID_DIM.z + (unsigned long)p.z;
		lex.fiber = f;
		destList->push_back(lex);

		p = p + v;
	}
}

void RasterizeFiberLex(rtsFiberNetwork* sourceNet, int f, list<LexID>* destList)
{
	/*resamples a fiber f and stores the resulting points as
	lexicographic values in LexicographicList.
	*/

	int num_points = sourceNet->FiberList[f].pointList.size();
	int p;

	point3D<float> pos = sourceNet->getNodeCoord(f, 0);

	point3D<float> coord;
	for(p=0; p<num_points; p++)
	{
		coord = sourceNet->getFiberPoint(f, p);
		RasterizeSegmentLex(f, pos, coord, destList);
		pos = coord;
	}
	RasterizeSegmentLex(f, pos, sourceNet->getNodeCoord(f, 1), destList);
}

list<LexID> RasterizeNetworkLex(rtsFiberNetwork* sourceNet, rtsDTGrid3D<float>* destGrid)
{
	/*This function creates a list of points on the focus network
	and stores them in a list as a lexicographic point:
	lp = z*sx*sy + y*sx + x
	*/

	//create a list that stores Lexicographic identifiers
	list<LexID> LexicographicList;

	int num_fibers = sourceNet->FiberList.size();
	int f;
	for(f=0; f<num_fibers; f++)
		RasterizeFiberLex(sourceNet, f, &LexicographicList);

	LexicographicList.sort();
	LexicographicList.unique();
	//LexicographicList.reverse();

	//nonLexVolume.SaveVOL("nonlexico.vol");

	//save the points in the focus field
	point3D<unsigned long> p;
	list<LexID>::iterator i;
	unsigned long lexID;
	for(i=LexicographicList.begin(); i!=LexicographicList.end(); i++)
	{
		lexID = (*i).ID;
		p.x = (lexID)/(GRID_DIM.y*GRID_DIM.z);
		p.y = (lexID - p.x*GRID_DIM.y*GRID_DIM.z)/GRID_DIM.z;
		p.z = (lexID - p.x*GRID_DIM.y*GRID_DIM.z - p.y*GRID_DIM.z);
		//cout<<"Pushing ";
		//p.print();
		//focusField->set(p.x, p.y, p.z, 255);
		destGrid->push(p.x, p.y, p.z, 255);
		//p.print();

	}
	return LexicographicList;
}

void EvaluateDistanceField(rtsDTGrid3D<float>* initialField)
{
	
	initialField->background = 255;
	(*initialField) = 0;
	FastSweeping3D(initialField, DILATION, VOXEL_SIZE);

}

float Gaussian(float x, float std)
{
	float prefix = 1.0/sqrt(2.0*3.14159*std*std);
	float suffix = exp(-(x*x)/(2.0*std*std));
	return prefix*suffix;
}

void EvaluateGaussianEnvelope(rtsDTGrid3D<float>* distanceField)
{
	
	rtsDTGrid3D<float>::iterator i;
	float G;
	float zero_value = Gaussian(0.0, STD);
	for(i=distanceField->begin(); i!=distanceField->after(); i++)
	{
		G = Gaussian(i.Value(), STD);
		//normalize the gaussian
		G /= zero_value;
		//cout<<G<<endl;
		i.SetValue(G);
	}
	distanceField->background = 0.0;

}

rtsDTGrid3D<float> CalculateDifference(rtsDTGrid3D<float>* grid0, rtsDTGrid3D<float>* grid1)
{
	rtsDTGrid3D<float> result;
	result = (*grid0) - (*grid1);

	
	rts_itkVolume<float> testVolume;
	testVolume.InsertDTGrid(&result);
	testVolume.SaveVOL("rasterized.vol");
	cout<<"Size: "<<testVolume.DimX()<<","<<testVolume.DimY()<<","<<testVolume.DimZ()<<endl;
	cout<<"Raw Size: "<<testVolume.DimX()*testVolume.DimY()*testVolume.DimZ()<<endl;
	

	return result;
}


float ComparePoints(float gold_val, float test_val, float std_1)
{
	/*if(gold_val > std_1)
		test_val = 0.0;
	if(test_val > std_1)
		gold_val = 0.0;
	*/
	if(test_val > std_1 && gold_val > std_1)
		return 0.0;
	return gold_val - test_val;

}

rtsDTGrid3D<float> CompareGrids(rtsDTGrid3D<float>* goldGrid, rtsDTGrid3D<float>* testGrid)
{
	rtsDTGrid3D<float> result;
	//create an iterator for each DT Grid
	rtsDTGrid3D<float>::iterator gold_i = goldGrid->begin();
	rtsDTGrid3D<float>::iterator test_i = testGrid->begin();

	//compute the value at one standard deviation
	float std_1 = Gaussian(STD, STD);
	cout<<"Value at 1 STD: "<<std_1<<endl;
	//create a variable to store the result value
	float result_value;

	//iterate both until one iterator has hit after()
	while(gold_i != goldGrid->after() && test_i != testGrid->after())
	{
		//if the iterators are at the same coordinate
		if(gold_i.Coord() == test_i.Coord())
		{
			//insert result of the comparison into the new grid
			result_value = ComparePoints(gold_i.Value(), test_i.Value(), std_1);
			result.push(gold_i.X1(), gold_i.X2(), gold_i.X3(), result_value);
			//increment both
			gold_i++; test_i++;
		}
		//add the lowest (lexicographically) value to the background, insert the result, and increment
		else if( (gold_i.Coord() < test_i.Coord()) )
		{
			result_value = ComparePoints(gold_i.Value(), testGrid->background, std_1);
			result.push(gold_i.X1(), gold_i.X2(), gold_i.X3(), result_value);
			gold_i++;
		}
		else if( (test_i.Coord() < gold_i.Coord()) )
		{
			result_value = ComparePoints(goldGrid->background, test_i.Value(), std_1);
			result.push(test_i.X1(), test_i.X2(), test_i.X3(), result_value);
			test_i++;
		}
	}

	cout<<"here"<<endl;

	//if the left iterator hasn't finished, iterate to finish it off
	while(gold_i != goldGrid->after())
	{
		result_value = ComparePoints(gold_i.Value(), testGrid->background, std_1);
		result.push(gold_i.X1(), gold_i.X2(), gold_i.X3(), result_value);
		gold_i++;
	}

	while(test_i != testGrid->after())
	{
		result_value = ComparePoints(goldGrid->background, test_i.Value(), std_1);
		result.push(test_i.X1(), test_i.X2(), test_i.X3(), result_value);
		test_i++;
	}


	return result;

	
	rts_itkVolume<float> testVolume;
	testVolume.InsertDTGrid(&result);
	testVolume.SaveVOL("rasterized.vol");
	cout<<"Size: "<<testVolume.DimX()<<","<<testVolume.DimY()<<","<<testVolume.DimZ()<<endl;
	cout<<"Raw Size: "<<testVolume.DimX()*testVolume.DimY()*testVolume.DimZ()<<endl;
	

	return result;
}

void IntegrateResult(rtsDTGrid3D<float>* inputGrid, float &total_positive, float &total_negative)
{
	total_positive = total_negative = 0.0;
	rtsDTGrid3D<float>::iterator i;
	float value;
	for(i = inputGrid->begin(); i != inputGrid->after(); i++)
	{
		value = i.Value();
		if(value > 0.0)
			total_positive += value;
		if(value < 0.0)
			total_negative += -value;
	}
}

float SumGrid(rtsDTGrid3D<float>* inputGrid)
{
	rtsDTGrid3D<float>::iterator i;
	float result = 0.0;
	for(i = inputGrid->begin(); i != inputGrid->after(); i++)
	{
		result += i.Value();
	}
	return result;
}

void StoreInFiber(rtsFiberNetwork* network, int fiber, float value)
{
	if(network->FiberList[fiber].FiberData == NULL)
	{
		network->FiberList[fiber].FiberData = (void*)(new float[2]);
		((float*)network->FiberList[fiber].FiberData)[0] = 0.0;
		((float*)network->FiberList[fiber].FiberData)[1] = 0.0;
	}

	((float*)network->FiberList[fiber].FiberData)[0] += value;
	((float*)network->FiberList[fiber].FiberData)[1] += 1.0;

}
void MapToExplicit(int radius, rtsDTGrid3D<float>* inGrid, list<LexID>* inLexList, rtsFiberNetwork* inNetwork)
{
	//create a template iterator to move over the input grid
	rtsMultiDirectionStencil3D<float> i;

	//create an iterator to move through the lexicographic list
	list<LexID>::iterator L = inLexList->begin();
	
	//set the stencil positions
	int x, y, z;
	for(x=-radius; x<=radius; x++)
		for(y=-radius; y<=radius; y++)
			for(z=-radius; z<=radius; z++)
			{
				i.addPosition(x, y, z);
			}

	i.Initialize(inGrid);
	unsigned long lex_id;
	float positives;
	int f;
	int p;
	for(i.StartPPP(); L != inLexList->end(); i.PPP())
	{
		//convert the iterator position to a Lexicographic ID
		lex_id = (unsigned long)i.X1()*GRID_DIM.y*GRID_DIM.z + (unsigned long)i.X2()*GRID_DIM.z + (unsigned long)i.X3();
		if(lex_id == (*L).ID)
		{
			//get the fiber ID associated with the point
			f = (*L).fiber;

			//get the positive value
			positives = 0.0;
			for(p = 0; p<i.numValues(); p++)
			{
				if(i.getValue(p) > 0)
					positives += i.getValue(p);
			}
			//if(positives > 0)
			//	cout<<positives<<endl;

			//store the value in the network
			StoreInFiber(inNetwork, f, positives);
			L++;

		}


	}
	


}
