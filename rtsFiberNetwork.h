#include <string>
#include <fstream>
#include <algorithm>
//#include "PerformanceData.h"
#include "objJedi.h"
#include "rtsPoint3d.h"
#include "rtsFilename.h"
#include <ANN/ANN.h>
//#include <exception>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
using namespace boost;


//Performance
//PerformanceData PD;

#ifndef RTS_FIBER_NETWORK_H
#define RTS_FIBER_NETWORK_H

//definitions for topologyNode labels
#define RTS_TOPOLOGY_NODE_NOEXIST	0
#define RTS_TOPOLOGY_NODE_INVALID	1
#define RTS_TOPOLOGY_NODE_VALID		2

#define RTS_TOPOLOGY_EDGE_NOEXIST	0
#define RTS_TOPOLOGY_EDGE_EXIST		1

//Properties for the topology graph
struct vertex_position_t
{
	typedef vertex_property_tag kind;
};
typedef property<vertex_position_t, point3D<float>, property<vertex_color_t, int> > VertexProperties;
typedef list<int> EdgeSequence;
typedef pair<EdgeSequence, EdgeSequence> EdgeMapping;
typedef vector<EdgeMapping> CoreGraphList;
typedef property<edge_weight_t, float, property<edge_color_t, EdgeSequence> > EdgeProperties;
typedef adjacency_list<listS, vecS, undirectedS, VertexProperties, EdgeProperties> TopologyGraph;
typedef graph_traits<TopologyGraph>::edge_descriptor TopologyEdge;
typedef graph_traits<TopologyGraph>::vertex_descriptor TopologyVertex;


EdgeSequence global_EdgeSequence;
float global_Weight;
list<float> global_NeighborWeights;
list<TopologyEdge> global_EdgeDescriptorSequence;
vector<TopologyVertex> global_Predecessors;
TopologyVertex global_Source;
pair<TopologyEdge, bool> BOOST_SmallestEdge(int v0, int v1, const TopologyGraph& G)
{
	pair<TopologyEdge, bool> edge_pair = edge(v0, v1, G);
	//if the edge doesn't exist, return the false pair
	if(!edge_pair.second)
		return edge_pair;

	//cout<<"Smallest edge: "<<endl;
	//BOOST_PrintGraph(G);

	//otherwise, make sure it is the edge with the least weight
	graph_traits<TopologyGraph>::out_edge_iterator oi, oi_end;
	TopologyEdge min_e = edge_pair.first;
	float min_weight = get(edge_weight_t(), G, min_e);

	for(tie(oi, oi_end) = out_edges(v0, G); oi!=oi_end; oi++)
	{
		if(target(*oi, G) == v1)
		{
			//cout<<"Edge weight for "<<v0<<": "<<get(edge_weight_t(), G, *oi)<<endl;
			if(get(edge_weight_t(), G, *oi) < min_weight)
			{
				min_weight = get(edge_weight_t(), G, *oi);
				min_e = *oi;
			}
		}
	}
	edge_pair.first = min_e;
	edge_pair.second = true;
	return edge_pair;
}
class boundary_bfs_visitor : public default_bfs_visitor
{
public:
	void discover_vertex(TopologyVertex v, const TopologyGraph& G) const throw(TopologyVertex)
	{
		//cout<<"discovered: "<<v<<endl;
		//global_EdgeSequence.push_back(get(vertex_color_t(), G, v));
		if(get(vertex_color_t(), G, v) >= 0 && v != global_Source)
		{
			TopologyVertex v0, v1;
			v0 = v1 = v;
			pair<TopologyEdge, bool> e;
			global_Weight = 0.0;
			while(v1 != global_Source)
			{
				v1 = global_Predecessors[v0];
				e = BOOST_SmallestEdge(v0, v1, G);
				global_EdgeSequence.push_back(get(edge_color_t(), G, e.first).front());
				global_EdgeDescriptorSequence.push_back(e.first);
				global_Weight += get(edge_weight_t(), G, e.first);
				v0 = v1;
			}
			//cout<<"Weight test: "<<global_Weight<<endl;
			throw v;
		}

	}
	void tree_edge(TopologyEdge e, const TopologyGraph& G) const
	{
		global_Predecessors[target(e, G)] = source(e, G);
	}
};



bool compare_edges(pair<TopologyEdge, float> e0, pair<TopologyEdge, float> e1)
{
	if(e0.second < e1.second)
		return true;
	else return false;
}


struct Fiber
{
	vector<point3D<float> > pointList;
	vector<float> errorList;
	int n0;
	int n1;
	float length;
	float error;
	int mapped_to;		//edge in another graph that maps to this one (-1 if a mapping doesn't exist)
	void* FiberData;

	Fiber()
	{
		error = 1.0;
		FiberData = NULL;
	}
};

struct Node
{
	point3D<float> p;
	void* NodeData;
	float error;
	int color;
	int incident;

	Node()
	{
		error = 1.0;
		NodeData = NULL;
		incident = 0;
	}
};

struct geometryPoint
{
	point3D<float> p;
	unsigned int fiberIdx;
	unsigned int pointIdx;
	float dist;
	int fiberID;
};

struct topologyEdge
{
	int label;
	unsigned int n0;
	unsigned int n1;
	float error;
};
struct topologyNode
{
	point3D<float> p;
	int label;
	list<unsigned int> connections;
	unsigned int compatible;
};




#define MAX_DIST	9999.0;

class rtsFiberNetwork
{
private:
	bool fiber_started;
	unsigned int num_points;
	float cull_value;	//used to cull fibers based on geometric error

	vector<geometryPoint> getNetPointSamples(float subdiv)
	{
		vector<geometryPoint> result;

		unsigned int f;
		list<point3D<float> > fiberPoints;
		list<point3D<float> >::iterator p;
		for(f = 0; f<FiberList.size(); f++)
		{
			fiberPoints.clear();
			fiberPoints = SubdivideFiber(f, subdiv);

			for(p=fiberPoints.begin(); p!=fiberPoints.end(); p++)
			{
				geometryPoint fp;
				fp.p = (*p);
				fp.fiberID = f;
				fp.dist = MAX_DIST;
				result.push_back(fp);
			}
		}
		return result;


	}

	void BF_ComputeL1Distance(vector<geometryPoint>* N0, vector<geometryPoint>* N1);
	void KD_ComputeEnvelopeDistance(rtsFiberNetwork* network, vector<geometryPoint>* Samples, float sigma)
	{
		//build the point arrays
		ANNpointArray dataPts0 = annAllocPts(Samples->size(), 3);
		for(unsigned int i=0; i<Samples->size(); i++)
		{
			dataPts0[i][0] = (*Samples)[i].p.x;
			dataPts0[i][1] = (*Samples)[i].p.y;
			dataPts0[i][2] = (*Samples)[i].p.z;
		}

		//create ANN variables
		ANNkd_tree* kdTree;
		ANNpoint queryPt = annAllocPt(3);
		ANNidxArray nearestIdx = new ANNidx[1];
		ANNdistArray nearestDist = new ANNdist[1];

		//compare network 0 to network 1
		//PD.StartTimer(LOG_N_DIST_BUILD0);
		kdTree = new ANNkd_tree(dataPts0, Samples->size(), 3);
		//PD.EndTimer(LOG_N_DIST_BUILD0);
		//PD.StartTimer(LOG_N_DIST_SEARCH0);

		//test each point in the network to the Samples list
		unsigned int f, p;
		int nodenum;
		float gauss_dist;
		for(f=0; f<network->FiberList.size(); f++)
		{
			//clear the error list
			network->FiberList[f].errorList.clear();

			//compute the distance at the nodes
			nodenum = network->FiberList[f].n0;
			queryPt[0] = network->NodeList[nodenum].p.x;
			queryPt[1] = network->NodeList[nodenum].p.y;
			queryPt[2] = network->NodeList[nodenum].p.z;
			kdTree->annkSearch(queryPt, 1, nearestIdx, nearestDist);
			gauss_dist = 1.0f - GaussianEnvelope(sqrtf((float)nearestDist[0]), sigma);
			network->NodeList[nodenum].error = gauss_dist;

			nodenum = network->FiberList[f].n1;
			queryPt[0] = network->NodeList[nodenum].p.x;
			queryPt[1] = network->NodeList[nodenum].p.y;
			queryPt[2] = network->NodeList[nodenum].p.z;
			kdTree->annkSearch(queryPt, 1, nearestIdx, nearestDist);
			gauss_dist = 1.0f - GaussianEnvelope(sqrtf((float)nearestDist[0]), sigma);
			network->NodeList[nodenum].error = gauss_dist;

			//compute the distance at each point along the fiber
			for(p=0; p<network->FiberList[f].pointList.size(); p++)
			{
				queryPt[0] = network->FiberList[f].pointList[p].x;
				queryPt[1] = network->FiberList[f].pointList[p].y;
				queryPt[2] = network->FiberList[f].pointList[p].z;
				kdTree->annkSearch(queryPt, 1, nearestIdx, nearestDist);
				gauss_dist = 1.0f - GaussianEnvelope(sqrtf((float)nearestDist[0]), sigma);
				network->FiberList[f].errorList.push_back(gauss_dist);
			}
		}
		/*for(int i=0; i<N1->size(); i++)
		{
			queryPt[0] = (*N1)[i].p.x;
			queryPt[1] = (*N1)[i].p.y;
			queryPt[2] = (*N1)[i].p.z;
			kdTree->annkSearch(queryPt, 1, nearestIdx, nearestDist);
			(*N1)[i].dist = sqrt(nearestDist[0]);
		}	*/
		//delete kdTree;
		//PD.EndTimer(LOG_N_DIST_SEARCH0);

		/*
		//compare network 1 to network 0
		PD.StartTimer(LOG_N_DIST_BUILD1);
		kdTree = new ANNkd_tree(dataPts1, N1->size(), 3);
		PD.EndTimer(LOG_N_DIST_BUILD1);
		PD.StartTimer(LOG_N_DIST_SEARCH1);
		for(int i=0; i<N0->size(); i++)
		{
			queryPt[0] = (*N0)[i].p.x;
			queryPt[1] = (*N0)[i].p.y;
			queryPt[2] = (*N0)[i].p.z;
			kdTree->annkSearch(queryPt, 1, nearestIdx, nearestDist);
			(*N0)[i].dist = sqrt(nearestDist[0]);
		}
		PD.EndTimer(LOG_N_DIST_SEARCH1);
		//delete kdTree;
		*/
		annClose();
		//delete kdTree;
	}

	void BD_ComputeL1Distance(vector<geometryPoint>* N0, vector<geometryPoint>* N1);
	list<point3D<float> > SubdivideSegment(point3D<float> p0, point3D<float> p1, float spacing)
	{

		//find the direction of travel
		vector3D<float> v = p1 - p0;

		//set the step size to the voxel size
		float length = v.Length();
		v.Normalize();

		float l;
		list<point3D<float> > result;
		point3D<float> p;
		for(l=0.0; l<length; l+=spacing)
		{
			p = p0 + v*l;
			result.push_back(p);
		}
		return result;


	}

	list<point3D<float> > SubdivideFiber(unsigned int f, float spacing)
	{
		list<point3D<float> > result;
		list<point3D<float> > segment;

		point3D<float> p0;
		point3D<float> p1;

		int node = FiberList[f].n0;
		p0 = NodeList[node].p;

		for(unsigned int p=0; p<FiberList[f].pointList.size(); p++)
		{
			segment.clear();
			//p0 = getFiberPoint(f, p-1);
			//p1 = getFiberPoint(f, p);
			p1 = FiberList[f].pointList[p];
			segment = SubdivideSegment(p0, p1, spacing);

			//result.push_back(p0);
			result.insert(result.end(), segment.begin(), segment.end());
			//result.push_back(p1);
			p0 = p1;
		}

		//subdivide the last segment
		node = FiberList[f].n1;
		p1 = NodeList[node].p;
		segment.clear();
		segment = SubdivideSegment(p0, p1, spacing);
		result.insert(result.end(), segment.begin(), segment.end());


		return result;
	}

	float GaussianEnvelope(float x, float std){return exp(-(x*x)/(2*std*std));}
	//vector<topologyNode> initNodeList(rtsFiberNetwork* network);
	void initTopologyGraph(vector<topologyNode>* Nodes, vector<topologyEdge>* Edges, rtsFiberNetwork* Network);
	void MapDeviationToNetwork(vector<geometryPoint>* source);
	void topLabelNodes(vector<topologyNode>* N0, vector<topologyNode>* N1, float epsilon);
	bool topMergeNode(vector<topologyNode>* NodeList, vector<topologyEdge>* EdgeList, unsigned int node);
	int topCollapse(vector<topologyNode>* Nodelist, vector<topologyEdge>* EdgeList);
	void topComputeMergePoints(vector<topologyNode>* NodeList);
	bool topDetectEdge(vector<topologyNode>* NodeList, vector<topologyEdge>* EdgeList, unsigned int node0, unsigned int node1);
	bool topDeleteEdge(vector<topologyNode>* NodeList, vector<topologyEdge>* EdgeList, unsigned int node0, unsigned int node1);
	unsigned int topGetNearestConnected(vector<topologyNode>* NodeList, unsigned int node, bool must_be_compatible);
	void ComputeBoundingVolume();
	float GeometryMetric(rtsFiberNetwork* network, float std)
	{
		//At this point each vertex in the network has an error in the range of [0 1]
		//This function computes the average error of each fiber and the entire network based
		//on the error at each vertex.
		unsigned int f, p;
		float fiber_length;
		float fiber_metric;
		float total_metric = 0.0;
		float total_length = 0.0;
		point3D<float> p0, p1;
		float e0, e1;
		int node;

		//compute the metric for every fiber in the network
		for(f=0; f<network->FiberList.size(); f++)
		{
			fiber_metric = 0.0;
			fiber_length = 0.0;

			//start at the first node
			node = network->FiberList[f].n0;
			p0 = network->NodeList[node].p;
			e0 = network->NodeList[node].error;

			//iterate through every intermediary node
			for(p=0; p<network->FiberList[f].errorList.size(); p++)
			{
				p1 = network->FiberList[f].pointList[p];
				e1 = network->FiberList[f].errorList[p];
				//keep a running sum of both the fiber length and the average error
				fiber_length += (p0 - p1).Length();
				fiber_metric += (e0+e1)/2 * (p0 - p1).Length();

				//shift
				p0 = p1;
				e0 = e1;

			}

			//end at the last fiber node
			node = network->FiberList[f].n1;
			p1 = network->NodeList[node].p;
			e1 = network->NodeList[node].error;
			fiber_length += (p0 - p1).Length();
			fiber_metric += (e0+e1)/2 * (p0 - p1).Length();

			//compute and store the average fiber error
			network->FiberList[f].error = fiber_metric/fiber_length;

			//keep a running total of the network error
			//do not include if the fiber is culled
			if(network->FiberList[f].error <= network->cull_value)
			{
				total_length += fiber_length;
				total_metric += fiber_metric;
			}
		}

		//compute the final error for the entire network
		return total_metric / total_length;
	}

	void MY_ComputeTopology(rtsFiberNetwork* testNetwork, float std);
	void BOOST_MapEdgesToFibers(rtsFiberNetwork* network, TopologyGraph& G)
	{
		/*
		//initialize the fiber mapping to -1
		for(int f=0; f<network->FiberList.size(); f++)
			network->FiberList[f].mapped_to = -1;

		//for each edge in the graph, set the mapping
		graph_traits<TopologyGraph>::edge_iterator ei, ei_end;
		int map;
		for(tie(ei, ei_end) = edges(G); ei != ei_end; ei++)
		{
			map = get(edge_color_t(), G, *ei);
			network->FiberList[map].mapped_to = 1;
		}
		*/

	}
	void BOOST_InvalidateValence2Vertices(TopologyGraph& G)
	{
		graph_traits<TopologyGraph>::vertex_iterator vi, vi_end;
		for(tie(vi, vi_end) = vertices(G); vi!= vi_end; vi++)
		{
			if(degree(*vi, G) == 2)
			{
				put(vertex_color_t(), G, *vi, -1);
				cout<<"invalidated"<<endl;
			}

		}
	}
	CoreGraphList BOOST_ComputeTopology(rtsFiberNetwork* testNetwork, float sigma)
	{
		//construct a graph representation of each network
		TopologyGraph GT = BOOST_BuildGraph(this);
		TopologyGraph T = BOOST_BuildGraph(testNetwork);

		cout<<"BEFORE CLEANING"<<endl;
		cout<<"GT****************************"<<endl;
		//BOOST_PrintGraph(GT);
		cout<<"T*****************************"<<endl;
		//BOOST_PrintGraph(T);

		//clean up degeneracies in the graphs
		BOOST_CleanGraph(GT, sigma);
		BOOST_CleanGraph(T, sigma);

		//color the vertices
		BOOST_KD_ColorGraph(GT, T, sigma);
		BOOST_KD_ColorGraph(T, GT, sigma);

		//color the vertices in GT based on the GT indices
		//the colors now match those in T
		BOOST_SaveGTIndices(GT);

		//apply the colors to the network
		BOOST_ApplyColorsToNetwork(GT, this);
		BOOST_ApplyColorsToNetwork(T, testNetwork);

		int false_positives = 0;
		int false_negatives = 0;

		//merge spurs
		//BOOST_PrintGraph(T);
		int false_positive_spines = BOOST_RemoveSpines(T);
		int false_negative_spines =  BOOST_RemoveSpines(GT);
		cout<<"False Positive Spines: "<<false_positive_spines<<endl;
		cout<<"False Negative Spines: "<<false_negative_spines<<endl;
		//BOOST_PrintGraph(T);


		false_positives += false_positive_spines;
		false_negatives += false_negative_spines;

		//merge invalid vertices together
		//this function also merges invalid vertices with valid vertices if their degree = 2

		//BOOST_PrintGraph(T);
		BOOST_MergeInvalidVertices(GT);
		BOOST_MergeInvalidVertices(T);
		//BOOST_PrintGraph(T);

		//merge the compatible vertices
		int compatible_spines = BOOST_MergeCompatibleVertices(T);
		cout<<"Compatibility errors in T: "<<compatible_spines<<endl;
		false_positives += compatible_spines;


		CoreGraphList core = BOOST_FindCore(GT, T);

		//BOOST_PrintGraph(GT);
		//BOOST_InvalidateValence2Vertices(GT);
		//BOOST_InvalidateValence2Vertices(T);
		//BOOST_PrintGraph(GT);
		//BOOST_MergeInvalidVertices(GT);
		//BOOST_MergeInvalidVertices(T);

		//core = BOOST_FindCore(GT, T);




		int T_topology_errors = num_edges(T) - core.size();
		int GT_topology_errors = num_edges(T) - core.size();

		false_positives += T_topology_errors;
		false_negatives += GT_topology_errors;

		cout<<"False Positive Edges: "<<false_positives<<endl;
		cout<<"False Negative Edges: "<<false_negatives<<endl;

		return core;
	}
	list<pair<TopologyVertex, EdgeSequence> > BOOST_FindNeighbors(TopologyGraph G, TopologyVertex node)
	{
		//Finds all colored vertices that can be reached from "node"
		pair<TopologyVertex, EdgeSequence> edge_map;
		list<pair<TopologyVertex, EdgeSequence> > result;

		do{
			global_Predecessors.clear();
			global_Predecessors.resize(num_vertices(G));
			global_Source = node;
			global_EdgeSequence.clear();
			global_EdgeDescriptorSequence.clear();
			global_Weight = 0.0;
			boundary_bfs_visitor vis;
			try{
				breadth_first_search(G, global_Source, visitor(vis));
			}
			catch(TopologyVertex& vert_id)
			{
				edge_map.first = vert_id;
				edge_map.second = global_EdgeSequence;
				result.push_back(edge_map);
				global_NeighborWeights.push_back(global_Weight);
				//clear_vertex(vert_id, G);
				remove_edge(global_EdgeDescriptorSequence.front(), G);
			}
		}while(!global_EdgeSequence.empty());
		return result;
	}

	TopologyGraph BOOST_FindCoreGraph(TopologyGraph G)
	{
		TopologyGraph result;
		//first insert all vertices
		graph_traits<TopologyGraph>::vertex_iterator vi, vi_end;
		graph_traits<TopologyGraph>::vertex_descriptor v;
		for(tie(vi, vi_end) = vertices(G); vi!=vi_end; vi++)
		{
			v = add_vertex(result);
			put(vertex_color_t(), result, v, get(vertex_color_t(), G, *vi));
			put(vertex_position_t(), result, v, get(vertex_position_t(), G, *vi));
		}

		//for each vertex in the graph, find all of the neighbors
		list<pair<TopologyVertex, EdgeSequence> > neighborhood;
		list<pair<TopologyVertex, EdgeSequence> >::iterator ni;
		list<float>::iterator wi;
		pair<TopologyEdge, bool> e;
		for(tie(vi, vi_end) = vertices(G); vi!=vi_end; vi++)
		{
			//only look for neighbors if the vertex has a color
			if(get(vertex_color_t(), G, *vi) >= 0)
			{
				global_NeighborWeights.clear();
				neighborhood = BOOST_FindNeighbors(G, *vi);
				for(ni = neighborhood.begin(), wi = global_NeighborWeights.begin(); ni != neighborhood.end(); ni++, wi++)
				{
					e = add_edge(*vi, (*ni).first, result);
					put(edge_color_t(), result, e.first, (*ni).second);
					//cout<<"Inserting weight: "<<*wi<<endl;
					put(edge_weight_t(), result, e.first, (*wi));

				}
				clear_vertex(*vi, G);
			}
		}
		//BOOST_PrintGraph(G);
		//BOOST_PrintGraph(result);
		return result;
	}

	CoreGraphList BOOST_CompareCoreGraphs(TopologyGraph GTc, TopologyGraph Tc)
	{
		//cout<<"Core Before:"<<endl;
		//BOOST_PrintGraph(GTc);

		CoreGraphList result;
		EdgeMapping TtoGT;
		pair<TopologyEdge, bool> e_Tc, e_GTc;
		graph_traits<TopologyGraph>::edge_iterator ei, ei_end, ei_temp;
		graph_traits<TopologyGraph>::vertex_descriptor v0, v1;

		//cout<<"Tc"<<endl;
		//BOOST_PrintGraph(Tc);
		//test each edge in T
		tie(ei, ei_end) = edges(Tc);
		while(ei!=ei_end)
		{
			//see if there is a corresponding edge in GT
			v0 = get(vertex_color_t(), Tc, source(*ei, Tc));
			v1 = get(vertex_color_t(), Tc, target(*ei, Tc));
			//e = edge(v0, v1, GTc);
			e_GTc = BOOST_SmallestEdge(v0, v1, GTc);
			//if(v0 == 1 || v1 == 1)
			//	cout<<"test"<<endl;

			if(e_GTc.second)
			{
				e_Tc = BOOST_SmallestEdge(source(*ei, Tc), target(*ei, Tc), Tc);

				//create the mapping
				TtoGT.first = get(edge_color_t(), Tc, e_Tc.first);
				TtoGT.second = get(edge_color_t(), GTc, e_GTc.first);
				result.push_back(TtoGT);

				//remove the edge
				remove_edge(e_GTc.first, GTc);
				remove_edge(e_Tc.first, Tc);
				tie(ei, ei_end) = edges(Tc);
			}
			else
				ei++;
		}



		//cout<<"Core after comparison: "<<endl;
		//cout<<"GTc"<<endl;
		//BOOST_PrintGraph(GTc);
		//cout<<"Tc"<<endl;
		//BOOST_PrintGraph(Tc);
		return result;
	}

	CoreGraphList NEW_ComputeTopology(rtsFiberNetwork* testNetwork, float sigma)
	{
		//cout<<"Number of Ground-Truth Nodes: "<<NodeList.size()<<endl;
		//cout<<"Number of Test-Case Nodes: "<<testNetwork->NodeList.size()<<endl;

		//construct a graph representation of each network
		TopologyGraph GT = BOOST_BuildGraph(this);
		TopologyGraph T = BOOST_BuildGraph(testNetwork);

		//color the vertices
		BOOST_KD_ColorGraphs(T, GT, sigma);

		//apply the colors to the network
		BOOST_ApplyColorsToNetwork(GT, this);
		BOOST_ApplyColorsToNetwork(T, testNetwork);

		//remove spurs
		//cout<<"Removing Spines...";
		//int false_positive_spines = BOOST_RemoveSpines(T);
		//int false_negative_spines =  BOOST_RemoveSpines(GT);
		//cout<<"done"<<endl;
		//cout<<"False Positive Spines: "<<false_positive_spines<<endl;
		//cout<<"False Negative Spines: "<<false_negative_spines<<endl;

		//initialize the global predecessor list
		//cout<<"Computing Core Connectivity for T...";
		TopologyGraph Tc = BOOST_FindCoreGraph(T);
		//cout<<"done"<<endl;
		//cout<<"Computing Core Connectivity for GT...";
		TopologyGraph GTc = BOOST_FindCoreGraph(GT);
		//cout<<"done"<<endl;
		//cout<<"Comparing Graphs...";
		return BOOST_CompareCoreGraphs(GTc, Tc);
		//cout<<"done"<<endl;
	}
	TopologyGraph BOOST_BuildGraph(rtsFiberNetwork* network)
	{
		//create the graph
		TopologyGraph g(network->NodeList.size());

		//for each fiber in the graph, create an edge
		int n0, n1;
		typedef std::pair<int, int> Edge;
		typedef std::pair<boost::graph_traits<TopologyGraph>::edge_descriptor, bool> EdgeIDType;
		EdgeIDType edge_id;
		Edge new_edge;
		EdgeProperties ep;
		EdgeSequence edge_sequence;
		for(unsigned int f=0; f<network->FiberList.size(); f++)
		{
			if(!network->isCulled(f))
			{
				n0 = network->FiberList[f].n0;
				n1 = network->FiberList[f].n1;
				new_edge = Edge(n0, n1);
				edge_id = add_edge(new_edge.first, new_edge.second, network->FiberList[f].error*network->FiberList[f].length, g);

				//add the starting edge color
				edge_sequence.clear();
				edge_sequence.push_back(f);
				put(edge_color_t(), g, edge_id.first, edge_sequence);
			}
			else
				cout<<"culled"<<endl;
		}

		//for each vertex in the graph, assign the position
		typedef property_map<TopologyGraph, vertex_position_t>::type PositionMap;
		PositionMap positions = get(vertex_position_t(), g);

		for(unsigned int v=0; v<network->NodeList.size(); v++)
			positions[v] = network->NodeList[v].p;

		return g;
	}
	void BOOST_KD_ColorGraph(TopologyGraph& G, TopologyGraph& compareTo, float sigma)
	{
		//get the number of vertices in compareTo
		int verts = num_vertices(compareTo);
		//allocate enough space in an ANN array
		ANNpointArray dataPts = annAllocPts(verts, 3);

		//get the vertex positions
		typedef property_map<TopologyGraph, vertex_position_t>::type PositionMap;
		PositionMap positions = get(vertex_position_t(), compareTo);
		//insert the positions into the ANN list
		for(int v=0; v<verts; v++)
		{
			dataPts[v][0] = positions[v].x;
			dataPts[v][1] = positions[v].y;
			dataPts[v][2] = positions[v].z;
		}
		//build the KD tree
		ANNkd_tree* kdTree;
		kdTree = new ANNkd_tree(dataPts, verts, 3);

		//PERFORM THE NEAREST NEIGHBOR SEARCH
		//allocate variables
		ANNpoint queryPt = annAllocPt(3);
		ANNidxArray nearestIdx = new ANNidx[1];
		ANNdistArray nearestDist = new ANNdist[1];

		//get the position map for G
		positions = get(vertex_position_t(), G);
		//get the vertex color map (valid/invalid)
		typedef property_map<TopologyGraph, vertex_color_t>::type ColorMap;
		ColorMap colors = get(vertex_color_t(), G);
		//get the vertex compatibility map
		//typedef property_map<TopologyGraph, vertex_compatibility_t>::type CompatibilityMap;
		//CompatibilityMap compatibility = get(vertex_compatibility_t(), G);
		//get the index property map
		//typedef property_map<TopologyGraph, vertex_index_t>::type IndexMap;
		//IndexMap index = get(vertex_index, G);

		//query each vertex in G
		typedef graph_traits<TopologyGraph>::vertex_iterator vertex_iter;
		std::pair<vertex_iter, vertex_iter> vp;
		point3D<float> pos;
		for (vp = vertices(G); vp.first != vp.second; ++vp.first)
		{
			pos = positions[*vp.first];
			queryPt[0] = pos.x;
			queryPt[1] = pos.y;
			queryPt[2] = pos.z;
			//perform the 1-NN search
			kdTree->annkSearch(queryPt, 1, nearestIdx, nearestDist);
			//if the distance is less than sigma, label as valid
			if(sqrt(nearestDist[0]) < sigma)
			{
				colors[*vp.first] = nearestIdx[0];
				//compatibility[*vp.first] = nearestIdx[0];
			}
			else
			{
				colors[*vp.first] = -1;
				//compatibility[*vp.first] = -1;
			}
		}
	}
	void BOOST_KD_ColorGraphs(TopologyGraph& T, TopologyGraph& GT, float sigma)
	{
		/*Colors both graphs for the connectivity metric:
		1) Create a kd-tree using the vertices in GT
		2) Find a vertex in GT near each vertex in T
			a) If a vertex in GT isn't found, the vertex in T is assigned a color of -1
			b) If a vertex in GT is found, the vertex in T is assigned the corresponding index
		3) Vertices in GT are assigned their own index if they are found (-1 if they aren't)
		*/

		//initialize each vertex in T with a color of -1
		graph_traits<TopologyGraph>::vertex_iterator vi, vi_end;
		for (tie(vi, vi_end) = vertices(T); vi != vi_end; vi++)
		{
			put(vertex_color_t(), T, *vi, -1);
		}

		//initialize each vertex in GT with -1
		for (tie(vi, vi_end) = vertices(GT); vi != vi_end; vi++)
		{
			put(vertex_color_t(), GT, *vi, -1);
		}
		//BOOST_PrintGraph(T);

		//CREATE THE KD TREE REPRESENTING GT
		//get the number of vertices in GT
		int verts = num_vertices(T);
		//allocate enough space in an ANN array
		ANNpointArray dataPts = annAllocPts(verts, 3);

		//get the vertex positions
		typedef property_map<TopologyGraph, vertex_position_t>::type PositionMap;
		PositionMap positions = get(vertex_position_t(), T);
		//insert the positions into the ANN list
		for(int v=0; v<verts; v++)
		{
			dataPts[v][0] = positions[v].x;
			dataPts[v][1] = positions[v].y;
			dataPts[v][2] = positions[v].z;
			//set the color for each vertex in GT to -1
			//put(vertex_color_t(), T, v, -1);
		}
		//build the KD tree
		ANNkd_tree* kdTree;
		kdTree = new ANNkd_tree(dataPts, verts, 3);

		//PERFORM THE NEAREST NEIGHBOR SEARCH
		//allocate variables
		ANNpoint queryPt = annAllocPt(3);
		ANNidxArray nearestIdx = new ANNidx[1];
		ANNdistArray nearestDist = new ANNdist[1];

		//get the position map for T
		//positions = get(vertex_position_t(), G);
		//get the vertex color map
		//typedef property_map<TopologyGraph, vertex_color_t>::type ColorMap;
		//ColorMap colors = get(vertex_color_t(), G);

		//query each vertex in T
		//std::pair<vertex_iter, vertex_iter> vp;
		point3D<float> pos;
		int num_undetected = 0;
		for (tie(vi, vi_end) = vertices(GT); vi != vi_end; vi++)
		{
			//pos = positions[*vp.first];
			pos = get(vertex_position_t(), GT, *vi);
			queryPt[0] = pos.x;
			queryPt[1] = pos.y;
			queryPt[2] = pos.z;
			//perform the 1-NN search
			kdTree->annkSearch(queryPt, 1, nearestIdx, nearestDist);
			//if the distance is less than sigma, label as valid
			if(sqrt(nearestDist[0]) < sigma)
			{
				//colors[*vp.first] = nearestIdx[0];
				//color T
				put(vertex_color_t(), T, nearestIdx[0], (int)*vi);
				//color GT
				put(vertex_color_t(), GT, *vi, (int)*vi);
				//compatibility[*vp.first] = nearestIdx[0];
			}
			else
			{
				//put(vertex_color_t(), T, *vi, -1);
				num_undetected++;
				//compatibility[*vp.first] = -1;
			}
		}

		//find undetected and falsely detected nodes
		//cout<<"Number of Undetected Nodes: "<<num_undetected<<endl;
		int false_positive = 0;
		int true_positive = 0;
		for (tie(vi, vi_end) = vertices(T); vi != vi_end; vi++)
		{
			if(get(vertex_color_t(), T, *vi) == -1)
				false_positive++;
			else
				true_positive++;
		}
		//cout<<"Number of False Positive Nodes: "<<false_positive<<endl;
		//cout<<"Number of True Positive Nodes: "<<true_positive<<endl;



	}
	void BOOST_CleanGraph(TopologyGraph& G, float sigma)
	{
		//This function cleans up degenerate cases in the ground truth and warns the user about them
		//These include:
		//		Ground truth nodes that are within sigma of other ground truth nodes

		//we do a 1-NN search across all vertices, looking for distance values < sigma

		bool error_detected = false;
		bool vertex_error;

		int merged_edges = 0;

		//go through each vertex
		graph_traits<TopologyGraph>::vertex_iterator vi, vtmp, vi_end;
		tie(vi, vi_end) = vertices(G);
		pair<graph_traits<TopologyGraph>::edge_descriptor, bool> e_remove;
		graph_traits<TopologyGraph>::adjacency_iterator ai, ai_end;
		graph_traits<TopologyGraph>::vertex_descriptor v_remove;
		point3D<float> vp, ap;
		vector3D<float> v_dist;
		float dist;
		while(vi != vi_end)
		{
			//get the color of the vertex
			vp = get(vertex_position_t(), G, *vi);

			//check every adjacent vertex
			vertex_error = false;
			for(tie(ai, ai_end) = adjacent_vertices(*vi, G); ai != ai_end; ai++)
			{
				//get the position of the adjacent vertex
				ap = get(vertex_position_t(), G, *ai);
				//find the distance between the two points
				v_dist = vp - ap;
				dist = v_dist.Length();

				if(dist <= sigma)
				{
					error_detected = true;
					vertex_error = true;
					v_remove = *vi;
					e_remove = edge(*vi, *ai, G);
				}
			}
			if(vertex_error)
			{
				//merge the vertices along edge "e_remove"
				BOOST_MergeVertices(G, e_remove.first, v_remove);

				//refresh the vertex iterator
				tie(vtmp, vi_end) = vertices(G);
			}
			else
				vi++;
		}


	}
	void BOOST_SaveGTIndices(TopologyGraph& G)
	{
		//changes the ground truth colors to vertex indices
		graph_traits<TopologyGraph>::vertex_iterator vi, vi_end;
		int index;
		for(tie(vi, vi_end) = vertices(G); vi != vi_end; vi++)
		{
			//only change the color if the vertex is valid
			if(get(vertex_color_t(), G, *vi) >= 0)
			{
				index = get(vertex_index_t(), G, *vi);
				put(vertex_color_t(), G, *vi, index);
			}
		}
	}
	void BOOST_ApplyColorsToNetwork(TopologyGraph& from, rtsFiberNetwork* to)
	{
		//get the index property
		typedef property_map<TopologyGraph, vertex_index_t>::type IndexMap;
		IndexMap index = get(vertex_index, from);
		//get the vertex color (valid/invalid)
		typedef property_map<TopologyGraph, vertex_color_t>::type ColorMap;
		ColorMap colors = get(vertex_color_t(), from);

		//go through each vertex and apply the appropriate color to the associated fiber node
		typedef graph_traits<TopologyGraph>::vertex_iterator vertex_iter;
		std::pair<vertex_iter, vertex_iter> vp;
		int idx;
		for (vp = vertices(from); vp.first != vp.second; ++vp.first)
		{
			idx = index[*vp.first];
			to->NodeList[idx].color = colors[*vp.first];
		}
		//cout << index[*vp.first] <<":"<<colors[*vp.first]<<"("<<compatibility[*vp.first]<<")= "<<positions[*vp.first]<<endl;

	}
	void BOOST_MergeVertices(TopologyGraph& G, graph_traits<TopologyGraph>::edge_descriptor e, graph_traits<TopologyGraph>::vertex_descriptor v)
	{
		//cout<<"MERGE VERTEX TEST-------------------------------"<<endl;
		//get the vertex descriptor for v0 and v1
		graph_traits<TopologyGraph>::vertex_descriptor v0, v1;
		v0 = v;
		v1 = target(e, G);
		if(v1 == v) v1 = source(e, G);

		//we will merge v0 to v1 (although the order doesn't matter)
		//well, okay, it does matter as far as the vertex position is concerned...maybe I'll fix that

		//get the current edge's colors
		EdgeSequence removed_colors = get(edge_color_t(), G, e);
		float removed_weight = get(edge_weight_t(), G, e);

		//first delete the current edge
		remove_edge(e, G);

		//for each edge coming into v0, create a corresponding edge into v1
		pair<graph_traits<TopologyGraph>::edge_descriptor, bool> new_edge;
		graph_traits<TopologyGraph>::in_edge_iterator ei, ei_end;
		EdgeSequence new_colors;
		float weight;
		for(tie(ei, ei_end) = in_edges(v0, G); ei != ei_end; ei++)
		{
			//create a new edge for each edge coming in to v0
			new_edge = add_edge(source(*ei, G), v1, G);
			new_colors = get(edge_color_t(), G, *ei);

			//add the removed colors to the new edge
			new_colors.insert(new_colors.end(), removed_colors.begin(), removed_colors.end());

			weight = get(edge_weight_t(), G, *ei) + removed_weight;
			put(edge_weight_t(), G, new_edge.first, weight);
			put(edge_color_t(), G, new_edge.first, new_colors);
		}

		//now remove v0 from the graph
		clear_vertex(v0, G);
		//remove_vertex(v0, G);

		//cout<<"END MERGE VERTEX TEST-------------------------------"<<endl;

	}
	void BOOST_MergeInvalidVertices(TopologyGraph& G)
	{
		//merges vertices that have color < 0

		//go through each vertex
		graph_traits<TopologyGraph>::vertex_iterator vi, vtmp, vi_end;
		tie(vi, vi_end) = vertices(G);
		int color;
		float min_weight;
		float weight;
		graph_traits<TopologyGraph>::edge_descriptor e_remove;
		graph_traits<TopologyGraph>::in_edge_iterator ei, ei_end;
		//graph_traits<TopologyGraph>::vertex_descriptor v_remove;
		bool vertex_removal;

		//merge neighbors that are both invalid
		//for each vertex
		for(tie(vi, vi_end) = vertices(G); vi!=vi_end; vi++)
		{
			vertex_removal = false;
			color = get(vertex_color_t(), G, *vi);

			//if the vertex is invalid
			if(color < 0)
			{

				//get the incident edge with the lowest weight
				min_weight = 99999.0;
				for(tie(ei, ei_end) = in_edges(*vi, G); ei != ei_end; ei++)
				{
					if(get(vertex_color_t(), G, source(*ei, G)) == get(vertex_color_t(), G, target(*ei, G)))
					{
						vertex_removal = true;
						weight = get(edge_weight_t(), G, *ei);
						cout<<"weight: "<<weight<<" source color: "<<color<<endl;
						if(weight <= min_weight)
						{
							min_weight = get(edge_weight_t(), G, *ei);
							e_remove = *ei;
						}
					}
				}
				//if a vertex is to be removed, remove it
				if(vertex_removal)
				{
					cout<<"Min Weight: "<<min_weight<<endl;
					//merge the vertices along edge "max_edge"
					BOOST_MergeVertices(G, e_remove, *vi);
				}

			}

		}

		//merge a vertex with its neighbor if it is invalid and degree=2
		for(tie(vi, vi_end) = vertices(G); vi!=vi_end; vi++)
		{
			//if the vertex is invalid
			if(color < 0 && degree(*vi, G) == 2)
			{
				tie(ei, ei_end) = in_edges(*vi, G);
				e_remove = *ei;
				BOOST_MergeVertices(G, e_remove, *vi);
			}
		}

	}
	int BOOST_RemoveSpines(TopologyGraph& G)
	{
		//merges vertices that have color < 0
		int merged_spines = 0;

		//go through each vertex
		graph_traits<TopologyGraph>::vertex_iterator vi, vtmp, vi_end;
		tie(vi, vi_end) = vertices(G);
		graph_traits<TopologyGraph>::edge_descriptor max_edge;
		graph_traits<TopologyGraph>::in_edge_iterator ei, ei_end;
		graph_traits<TopologyGraph>::vertex_descriptor v_next, v_this;
		graph_traits<TopologyGraph>::adjacency_iterator ai, ai_end;
		int v_degree;

		//for each vertex
		for(tie(vi, vi_end) = vertices(G); vi!=vi_end; vi++)
		{
			//if the vertex is invalid AND degree=1
			v_degree = degree(*vi, G);
			if(get(vertex_color_t(), G, *vi) < 0 && v_degree == 1)
			{
				//delete the entire spur
				v_this = *vi;
				do{
					//get the adjacent vertex
					tie(ai, ai_end) = adjacent_vertices(v_this, G);
					//cout<<*ai<<endl;
					v_next = *ai;
					//clear the current vertex (there should only be one edge)
					clear_vertex(v_this, G);
					v_this = v_next;
				}while(degree(v_this, G) == 1 && get(vertex_color_t(), G, v_this) < 0);
				merged_spines++;
			}
		}


		return merged_spines;
	}
	int BOOST_MergeCompatibleVertices(TopologyGraph& G)
	{
		int merged_edges = 0;

		//go through each vertex
		graph_traits<TopologyGraph>::vertex_iterator vi, vtmp, vi_end;
		tie(vi, vi_end) = vertices(G);
		int color;
		pair<graph_traits<TopologyGraph>::edge_descriptor, bool> e_remove;
		graph_traits<TopologyGraph>::adjacency_iterator ai, ai_end;
		graph_traits<TopologyGraph>::vertex_descriptor v_remove;
		int num_incident;
		bool merge_found;
		while(vi != vi_end)
		{
			//get the color of the vertex
			color = get(vertex_color_t(), G, *vi);
			//if the vertex is valid
			if(color >= 0)
			{
				//run through each adjacent vertex
				merge_found = false;
				for(tie(ai, ai_end) = adjacent_vertices(*vi, G); ai != ai_end; ai++)
				{
					//if the two vertices are compatible
					if(color == get(vertex_color_t(), G, *ai) && (*vi != *ai))
					{
						v_remove = *vi;
						e_remove = edge(*vi, *ai, G);
						merge_found = true;
					}
				}
				if(merge_found)
				{
					num_incident = degree(v_remove, G);
					BOOST_MergeVertices(G, e_remove.first, v_remove);
					if(num_incident == 1) merged_edges++;
					tie(vtmp, vi_end) = vertices(G);
				}
				else vi++;
			}
			else
				vi++;
		}

		return merged_edges;
	}

	EdgeSequence BOOST_RemoveEdge(int v0, int v1, TopologyGraph& G)
	{
		//look for a direct edge
		pair<graph_traits<TopologyGraph>::edge_descriptor, bool> e0, e1;
		e0 = BOOST_SmallestEdge(v0, v1, G);

		//if it exists, remove it and return the edge sequence
		EdgeSequence sequence0, sequence1;
		sequence0.clear();
		if(e0.second)
		{
			sequence0 = get(edge_color_t(), G, e0.first);
			remove_edge(e0.first, G);
		}
		//otherwise look for an indirect edge

		else
		{
			graph_traits<TopologyGraph>::adjacency_iterator ai, ai_end;
			//for each vertex adjacent to v0
			for(tie(ai, ai_end) = adjacent_vertices(v0, G); ai!=ai_end; ai++)
			{
				//if the adjacent vertex is invalid
				if(get(vertex_color_t(), G, *ai) < 0)
				{
					//see if v1 is also connected
					e1 = BOOST_SmallestEdge(v1, *ai, G);
					//if it is
					if(e1.second)
					{
						sequence0 = get(edge_color_t(), G, e1.first);
						e0 = BOOST_SmallestEdge(v0, *ai, G);
						sequence1 = get(edge_color_t(), G, e0.first);
						sequence0.insert(sequence0.end(), sequence1.begin(), sequence1.end());
						remove_edge(e0.first, G);
						remove_edge(e1.first, G);
						return sequence0;
					}
				}
			}
		}
		return sequence0;

	}
	CoreGraphList BOOST_FindCore(TopologyGraph& GT, TopologyGraph& T)
	{
		CoreGraphList Gc;

		int FP = 0;
		int FN = 0;

		//Find edges that exist in both GT and T
		EdgeSequence T_Path, GT_Path;
		EdgeMapping path_pair;

		//create iterators for the edges in T
		graph_traits<TopologyGraph>::edge_iterator ei, ei_end;

		//for each edge in T, remove the corresponding edge in GT (if it exists) and put the details in the core garph list
		int v0, v1;
		pair<graph_traits<TopologyGraph>::edge_descriptor, bool> e;
		pair<graph_traits<TopologyGraph>::edge_descriptor, bool> e_T, e_GT;

		graph_traits<TopologyGraph>::vertex_descriptor s, t;
		list<graph_traits<TopologyGraph>::edge_descriptor> to_remove;
		list<graph_traits<TopologyGraph>::edge_descriptor>::iterator remove_iter;

		/*//add all edges to a list
		list<pair<TopologyEdge, float>> EdgeList;
		for(tie(ei, ei_end) = edges(T); ei!=ei_end; ei++)
		{
			pair<TopologyEdge, float> edge_data;
			edge_data.first = *ei;
			edge_data.second = get(edge_weight_t(), T, *ei);
			EdgeList.push_back(edge_data);
		}
		EdgeList.sort(compare_edges);*/

		//find direct connections
		for(tie(ei, ei_end) = edges(T); ei!=ei_end; ei++)
		{
			//get the ID for the corresponding vertices in GT
			s = source(*ei, T);
			t = target(*ei, T);
			v0 = get(vertex_color_t(), T, s);
			v1 = get(vertex_color_t(), T, t);

			//if both vertices are valid
			if(v0 >= 0 && v1 >= 0)
			{
				GT_Path = BOOST_RemoveEdge(v0, v1, GT);
				//if there was a corresponding edge in GT
				if(GT_Path.size() > 0)
				{
					//get the edge path from T
					T_Path = get(edge_color_t(), T, *ei);
					//create the path pair and add it to the core graph list
					path_pair.first = T_Path;
					path_pair.second = GT_Path;
					Gc.push_back(path_pair);
				}
				else
					FP++;
				//mark the edge for removal from T
				to_remove.push_back(*ei);

			}
		}
		for(remove_iter = to_remove.begin(); remove_iter != to_remove.end(); remove_iter++)
			remove_edge(*remove_iter, T);


		//find indirect connections (these are connections that use one invalid point)
		graph_traits<TopologyGraph>::adjacency_iterator ai, ai_end;
		bool match_found;
		while(num_edges(T) > 0)
		{
			tie(ei, ei_end) = edges(T);

			//s = valid point, t = invalid
			s = source(*ei, T);
			if(get(vertex_color_t(), T, s) < 0)
			{
				t = s;
				s = target(*ei, T);
			}
			else
				t = target(*ei, T);
			//for each point adjacent to the invalid point
			match_found = false;
			v0 = get(vertex_color_t(), T, s);
			for(tie(ai, ai_end) = adjacent_vertices(t, T); ai!=ai_end; ai++)
			{
				//get the edge path from T
				e = edge(t, *ai, T);
				//make sure that the first and second edges are not the same
				if(e.first != *ei)
				{
					v1 = get(vertex_color_t(), T, *ai);
					GT_Path = BOOST_RemoveEdge(v0, v1, GT);
					//if there was a corresponding edge in GT
					if(GT_Path.size() > 0)
					{
						match_found = true;
						T_Path = get(edge_color_t(), T, e.first);
						EdgeSequence temp = get(edge_color_t(), T, *ei);
						T_Path.insert(T_Path.end(), temp.begin(), temp.end());

						//create the path pair and add it to the core graph list
						path_pair.first = T_Path;
						path_pair.second = GT_Path;
						Gc.push_back(path_pair);

						//remove both edges and break
						remove_edge(e.first, T);
						remove_edge(*ei, T);
						break;
					}
				}
			}
			if(!match_found)
			{
				FP++;
				remove_edge(*ei, T);
			}
		}

		cout<<"False Positive Edges in Core: "<<FP<<endl;



		return Gc;
	}

	TopologyGraph BOOST_RemapGraph(TopologyGraph& G, rtsFiberNetwork* network);
	void BOOST_PrintGraph(TopologyGraph G)
	{
		//get the index property
		typedef property_map<TopologyGraph, vertex_index_t>::type IndexMap;
		IndexMap index = get(vertex_index, G);

		//get the position property
		typedef property_map<TopologyGraph, vertex_position_t>::type PositionMap;
		PositionMap positions = get(vertex_position_t(), G);
		//get the vertex color (valid/invalid)
		typedef property_map<TopologyGraph, vertex_color_t>::type ColorMap;
		ColorMap colors = get(vertex_color_t(), G);
		//get vertex compatibility
		//typedef property_map<TopologyGraph, vertex_compatibility_t>::type CompatibilityMap;
		//CompatibilityMap compatibility = get(vertex_compatibility_t(), G);

		std::cout << "vertices(g) = "<<endl;
		typedef graph_traits<TopologyGraph>::vertex_iterator vertex_iter;
		std::pair<vertex_iter, vertex_iter> vp;
		for (vp = vertices(G); vp.first != vp.second; ++vp.first)
			cout << index[*vp.first] <<"("<<colors[*vp.first]<<") = "<<positions[*vp.first]<<endl;
		cout << endl;

		//get the edge weight property
		typedef property_map<TopologyGraph, edge_weight_t>::type WeightMap;
		WeightMap weights = get(edge_weight_t(), G);

		std::cout << "edges(source, dest): weight to_fiber"<<endl;
		graph_traits<TopologyGraph>::edge_iterator ei, ei_end;
		EdgeSequence edge_sequence;
		for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei)
		{
			std::cout << "(" << index[source(*ei, G)] << "," << index[target(*ei, G)] << "): "<<weights[*ei]<<" {";
			edge_sequence = get(edge_color_t(), G, *ei);
			for(EdgeSequence::iterator i = edge_sequence.begin(); i != edge_sequence.end(); i++)
				cout<<*i<<" ";
			cout<<"}"<<endl;
		}

		std::cout << std::endl;


	}

	//network geometry functions
	double computeFiberLength(int f)
	{
		//return the length of the fiber f
		point3D<float> p0 = NodeList[FiberList[f].n0].p;
		point3D<float> p1;

		double length = 0.0;
		int num_points = FiberList[f].pointList.size();
		int p;
		for(p=0; p<num_points; p++)
		{
			p1 = FiberList[f].pointList[p];
			length += (p1 - p0).Length();
			p0 = p1;
		}
		p1 = NodeList[FiberList[f].n1].p;
		length += (p1 - p0).Length();
		return length;

	}
	void refreshFiberLengths()
	{
		for(unsigned int f=0; f<FiberList.size(); f++)
			FiberList[f].length = (float)computeFiberLength(f);
	}

	void refreshIncidence()
	{
		for(unsigned int n=0; n<NodeList.size(); n++)
			NodeList[n].incident = 0;

		for(unsigned int f=0; f<FiberList.size(); f++)
		{
			NodeList[FiberList[f].n0].incident++;
			NodeList[FiberList[f].n1].incident++;
		}
	}

public:
	vector<Fiber> FiberList;
	vector<Node> NodeList;
	point3D<float> min_pos;
	point3D<float> max_pos;

	//constructors
	rtsFiberNetwork()
	{
		fiber_started = false;
		num_points = false;
		cull_value = 1.0;
		min_pos = point3D<float>(9999, 9999, 9999);
		max_pos = point3D<float>(-9999, -9999, -9999);
	}

	//get functions
	point3D<float> getNodeCoord(int node){return NodeList[node].p;}
	point3D<float> getNodeCoord(int fiber, bool node);
	point3D<float> getFiberPoint(unsigned int fiber, unsigned int point);

	//network culling
	float getCullValue(){return cull_value;}
	void setCullValue(float cull){cull_value = cull;}
	bool isCulled(int f)
	{
		if(FiberList[f].error <= cull_value)
			return false;
		else return true;
	}

	//statistics functions
	double getTotalLength();
	double getFiberLength(int f)
	{
		return FiberList[f].length;
	}



	//drawing functions
	void StartFiber(float x, float y, float z)
	{
		fiber_started = true;
		num_points = 1;

		//create a start node and an end node
		Node n;
		n.p = point3D<float>(x, y, z);
		NodeList.push_back(n);
		NodeList.push_back(n);

		//create a fiber
		Fiber f;
		f.n0 = NodeList.size()-2;
		f.n1 = NodeList.size()-1;
		FiberList.push_back(f);

	}
	void StartFiber(int node)
	{
		fiber_started = true;
		num_points = 1;

		//the start node is specified
		//specify the end node
		Node n = NodeList[node];
		NodeList.push_back(n);

		//create a fiber
		Fiber f;
		f.n0 = node;
		f.n1 = NodeList.size()-1;
		FiberList.push_back(f);
	}

	void ContinueFiber(float x, float y, float z)
	{
 		if(!fiber_started)
		{
			StartFiber(x, y, z);
			return;
		}
		num_points++;

		//store the last node coordinate in the fiber list
		int f = FiberList.size() - 1;
		int n = NodeList.size() - 1;

		if(num_points > 2)
			FiberList[f].pointList.push_back(NodeList[n].p);
		NodeList[n].p = point3D<float>(x, y, z);
		//NodeList[n].p.print();
		//cout<<endl;
	}

	void EndFiber()
	{
		fiber_started = false;
	}
	void StartBranch(int node)
	{

		num_points = 1;
		fiber_started = true;

		//create an end node only
		Fiber f;
		f.n0 = node;

		Node n = NodeList[node];
		NodeList.push_back(n);


		f.n1 = NodeList.size()-1;
		FiberList.push_back(f);

	}

	void EndBranch(int node)
	{
		int n = NodeList.size() - 1;
		int f = FiberList.size() - 1;

		//store the current node position in the fiber list
		FiberList[f].pointList.push_back(NodeList[n].p);

		//set the picked node as the fiber destination
		FiberList[f].n1 = node;

		//remove the final node (now not used)
		NodeList.pop_back();

		fiber_started = false;
	}

	void ConnectFiber(unsigned int node)
	{
		if(!fiber_started)
		{
			StartFiber(node);
			return;
		}
		if(node >= NodeList.size() - 1)
			return;

		int f = FiberList.size() - 1;

		//add the last point
		FiberList[f].pointList.push_back(NodeList[FiberList[f].n1].p);

		//attach to the specified node
		FiberList[f].n1 = node;
		NodeList.pop_back();
		fiber_started = false;
	}
	void StartFiber(point3D<float> p){StartFiber(p.x, p.y, p.z);}
	void ContinueFiber(point3D<float> p){ContinueFiber(p.x, p.y, p.z);}

	void Clear()
	{
		FiberList.clear();
		NodeList.clear();
		fiber_started = false;
	}

	//loading functions
	void MergeNodes(float sigma)
	{
		//create a KD tree for all nodes in the network
		//get the number of vertices in GT
		int verts = NodeList.size();
		//allocate enough space in an ANN array
		ANNpointArray dataPts = annAllocPts(verts, 3);

		//insert the positions into the ANN list
		for(int v=0; v<verts; v++)
		{
			dataPts[v][0] = NodeList[v].p.x;
			dataPts[v][1] = NodeList[v].p.y;
			dataPts[v][2] = NodeList[v].p.z;
		}
		//build the KD tree
		ANNkd_tree* kdTree;
		kdTree = new ANNkd_tree(dataPts, verts, 3);

		//create a map
		vector<int> NodeMap;
		NodeMap.resize(NodeList.size(), -1);

		ANNpoint queryPt = annAllocPt(3);
		ANNidxArray nearestIdx = new ANNidx[1];
		ANNdistArray nearestDist = new ANNdist[1];

		for(unsigned int n=0; n<NodeList.size(); n++)
		{
			queryPt[0] = NodeList[n].p.x;
			queryPt[1] = NodeList[n].p.y;
			queryPt[2] = NodeList[n].p.z;
			//perform the 1-NN search
			kdTree->annkSearch(queryPt, 1, nearestIdx, nearestDist);

			NodeMap[n] = nearestIdx[0];
		}

		//set the nodes for each fiber to those mapped
		for(unsigned int f=0; f<FiberList.size(); f++)
		{
			FiberList[f].n0 = NodeMap[FiberList[f].n0];
			FiberList[f].n1 = NodeMap[FiberList[f].n1];
		}

		RemoveExcessNodes();
		refreshIncidence();
	}
	void RemoveExcessNodes()
	{
		vector<int> NodeMap;
		NodeMap.resize(NodeList.size(), -1);
		vector<int> FiberNum;
		FiberNum.resize(NodeList.size(), 0);

		vector<Node> newNodeList;

		//run through each fiber
		int node;
		for(unsigned int f=0; f<FiberList.size(); f++)
		{
			node = FiberList[f].n0;
			//if this node has not been encountered
			if(NodeMap[node] == -1)
			{
				NodeMap[node] = newNodeList.size();
				FiberList[f].n0 = NodeMap[node];
				newNodeList.push_back(NodeList[node]);

				FiberNum[node]++;
			}
			else
			{
				FiberList[f].n0 = NodeMap[node];

				FiberNum[node]++;
			}

			node = FiberList[f].n1;
			if(NodeMap[node] == -1)
			{
				NodeMap[node] = newNodeList.size();
				FiberList[f].n1 = NodeMap[node];
				newNodeList.push_back(NodeList[node]);

				FiberNum[node]++;
			}
			else
			{
				FiberList[f].n1 = NodeMap[node];

				FiberNum[node]++;
			}
		}

		NodeList = newNodeList;
	}
	void LoadSWC(string filename)
	{
		//open the file
		ifstream infile;
		infile.open(filename.c_str());
		if(!infile.is_open())
			return;

		//variables to read
		int id, type, parent;
		int prev_id = 0;
		float x, y, z, r;

		vector<int> branchIDList;
		vector<int> correspondingNodeList;

		vector<int>::iterator iter;
		int index;

		char c;
		//first pass, get branch points and find bounds
		while(!infile.eof())
		{
			c = infile.peek();
			if((c <'0' || c > '9') && c != 32)
				infile.ignore(9999, '\n');
			else
			{
				infile>>id;
				infile>>type;
				infile>>x;
				infile>>y;
				infile>>z;
				infile>>r;
				infile>>parent;


				//root nodes and branch nodes
				if(parent != id-1 && parent != -1)
					branchIDList.push_back(parent);
				if(x < min_pos.x) min_pos.x = x;
				if(y < min_pos.y) min_pos.y = y;
				if(z < min_pos.z) min_pos.z = z;
				if(x > max_pos.x) max_pos.x = x;
				if(y > max_pos.y) max_pos.y = y;
				if(z > max_pos.z) max_pos.z = z;
			}
		}//end while
		//sort the branch points
		sort(branchIDList.begin(), branchIDList.end());
		//set the number of corresponding nodes
		correspondingNodeList.resize(branchIDList.size());

		//second pass, build the tree
		infile.close();
		infile.clear();
		infile.open(filename.c_str());
		while(!infile.eof())
		{
			c = infile.peek();
			if((c <'0' || c > '9') && c != 32)
				infile.ignore(9999, '\n');
			else
			{
				infile>>id;
				infile>>type;
				infile>>x;
				infile>>y;
				infile>>z;
				infile>>r;
				infile>>parent;


				//if we are starting a new fiber
				if(parent == -1)
					StartFiber(x, y, z);
				else
				{
					//see if the parent is in the branch list
					iter = find(branchIDList.begin(), branchIDList.end(), parent);
					//if the parent point is a branch point
					if(iter != branchIDList.end())
					{
						index = iter - branchIDList.begin();
						StartBranch(correspondingNodeList[index]);
						ContinueFiber(x, y, z);
					}
					else
						ContinueFiber(x, y, z);
				}
				//if the current node is in the branch list
				iter = find(branchIDList.begin(), branchIDList.end(), id);
				if(iter != branchIDList.end())
				{
					//get the id and store the corresponding node value
					index = iter - branchIDList.begin();
					correspondingNodeList[index] = NodeList.size()-1;
				}
			}
		}//end while
		refreshFiberLengths();
		refreshIncidence();

	}

	void LoadOBJ(string filename)
	{
		//first load the OBJ file
		rtsOBJ objFile;
		objFile.LoadFile(filename.c_str());

		/*Create two lists. Each element represents a point in the OBJ file.  We will
		first go through each fiber (line) and find the vertex associated with the two ends of the fiber.
		The validList value is "true" if the associated point in the OBJ file is a junction.  It is "false"
		if the point is just an intermediate point on a fiber. The pointList value stores the new value of the junction
		in the NodeList of this FiberNetwork structure.*/

		vector<bool> validList;
		validList.assign(objFile.v_list.size(), false);
		vector<unsigned int> pointList;
		pointList.assign(objFile.v_list.size(), 0);

		/*Run through each fiber:
		1) See if each fiber node has already been created by looking at validList (true = created)
		2) If the node has not been created, create it and set validList to true and pointList to the node index in this structure.
		3) Create an empty fiber with the appropriate node values assigned.*/
		unsigned int f;
		unsigned int line_verts;
		unsigned int v0, vn;
		Node node_new;
		for(f=0; f<objFile.getNumLines(); f++)
		{
			//find out how many vertices there are in the line
			line_verts = objFile.getNumLineVertices(f);
			v0 = objFile.getLineVertex(f, 0);
			vn = objFile.getLineVertex(f, line_verts - 1);

			//if the nodes don't exist, create them
			if(!validList[v0])
			{
				node_new.p = objFile.getVertex3d(v0);
				pointList[v0] = NodeList.size();
				NodeList.push_back(node_new);
				validList[v0] = true;
			}
			if(!validList[vn])
			{
				node_new.p = objFile.getVertex3d(vn);
				pointList[vn] = NodeList.size();
				NodeList.push_back(node_new);
				validList[vn] = true;
			}
			//create the new fiber
			Fiber fiber_new;
			fiber_new.n0 = pointList[v0];
			fiber_new.n1 = pointList[vn];

			//get all of the intermediate line vertices and insert them in the fiber
			for(unsigned int i=1; i<line_verts-1; i++)
			{
				fiber_new.pointList.push_back(objFile.getVertex3d(objFile.getLineVertex(f, i)));
				//fiber_new.pointList.push_back(objFile.get
			}
			FiberList.push_back(fiber_new);
		}

		RemoveExcessNodes();

		ComputeBoundingVolume();
		refreshFiberLengths();
		refreshIncidence();

	}
	void LoadFile(string filename)
	{
		rtsFilename file = filename;
		string extension = file.getExtension();
		if(extension.compare("obj") == 0)
			LoadOBJ(filename);
		else if(extension.compare("swc") == 0)
			LoadSWC(filename);

	}


	//saving functions
	void SaveOBJ(string filename)
	{
		ofstream outfile;
		outfile.open(filename.c_str());

		//output all vertices

		//first output all nodes
		for(unsigned int n=0; n<NodeList.size(); n++)
			outfile<<"v "<<NodeList[n].p.x<<" "<<NodeList[n].p.y<<" "<<NodeList[n].p.z<<endl;
		//then output all fiber points
		for(unsigned int f=0; f<FiberList.size(); f++)
			for(unsigned int p=0; p<FiberList[f].pointList.size(); p++)
				outfile<<"v "<<FiberList[f].pointList[p].x<<" "<<FiberList[f].pointList[p].y<<" "<<FiberList[f].pointList[p].z<<endl;

		//now output each of the fibers
		int i = NodeList.size() + 1;
		for(unsigned int f=0; f<FiberList.size(); f++)
		{
			outfile<<"l ";
			outfile<<FiberList[f].n0+1<<" ";
			for(unsigned int p=0; p<FiberList[f].pointList.size(); p++)
			{
				outfile<<i<<" ";
				i++;
			}
			outfile<<FiberList[f].n1+1;
			outfile<<endl;
		}

		outfile.close();

	}

	//transform functions
	void Translate(float x, float y, float z);
	void Translate(point3D<float> p){Translate(p.x, p.y, p.z);}

	void Oscillate(float frequency, float magnitude)
	{
		//impliment a sinusoidal oscillation along each fiber
		vector3D<float> side(0.0, 0.0, 1.0);
		vector3D<float> dir, normal;
		point3D<float> p0, p1;
		float t;

		//for each fiber
		for(unsigned int f=0; f<FiberList.size(); f++)
		{
			p0 = NodeList[FiberList[f].n0].p;
			t=0.0;

			num_points = FiberList[f].pointList.size();
			for(unsigned int p = 0; p<num_points; p++)
			{
				p1 = FiberList[f].pointList[p];
				dir = p0 - p1;
				t+= dir.Length();
				normal = dir.X(side);
				normal.Normalize();
				FiberList[f].pointList[p] = FiberList[f].pointList[p] + magnitude*sin(t*frequency)*normal;

				p0 = p1;
			}


		}

	}

	void Crop(float px, float py, float pz, float sx, float sy, float sz)
	{
		vector<Fiber> newFiberList;
		for(unsigned int f=0; f<FiberList.size(); f++)
		{
			int n0 = FiberList[f].n0;
			int n1 = FiberList[f].n1;
			point3D<float> p0 = NodeList[n0].p;
			point3D<float> p1 = NodeList[n1].p;

			if(p0.x > px && p0.y > py && p0.z > pz &&
			   p1.x > px && p1.y > py && p1.z > pz &&
			   p0.x < px+sx && p0.y < py+sy && p0.z < pz+sz &&
			   p1.x < px+sx && p1.y < py+sy && p1.z < pz+sz)
			{
				newFiberList.push_back(FiberList[f]);
			}
		}
		FiberList.clear();
		FiberList = newFiberList;
		RemoveExcessNodes();
		ComputeBoundingVolume();
	}
	void ThresholdLength(float min_length, float max_length = 99999)
	{
		vector<Fiber> newFiberList;
		refreshFiberLengths();
		for(unsigned int f=0; f<FiberList.size(); f++)
		{
			if(FiberList[f].length > min_length && FiberList[f].length < max_length)
			{
				newFiberList.push_back(FiberList[f]);
			}
		}
		FiberList.clear();
		FiberList = newFiberList;
		RemoveExcessNodes();
		ComputeBoundingVolume();
	}
	void ThresholdSpines(float min_length)
	{
		vector<Fiber> newFiberList;
		refreshIncidence();
		refreshFiberLengths();
		for(unsigned int f=0; f<FiberList.size(); f++)
		{
			if(FiberList[f].length > min_length || (NodeList[FiberList[f].n0].incident > 1 && NodeList[FiberList[f].n1].incident > 1))
			{
				newFiberList.push_back(FiberList[f]);
			}
		}
		FiberList.clear();
		FiberList = newFiberList;
		RemoveExcessNodes();
		ComputeBoundingVolume();

	}
	//subdivision
	void SubdivideNetwork(float spacing)
	{
		list<point3D<float> > subdivided;
		list<point3D<float> >::iterator p;
		for(unsigned int f=0; f<FiberList.size(); f++)
		{
			//get the subdivided fiber
			subdivided.clear();
			subdivided = SubdivideFiber(f, spacing);

			//clean up the current fiber
			FiberList[f].pointList.clear();
			//copy the subdivided fiber into the current fiber point list
			for(p = subdivided.begin(); p!=subdivided.end(); p++)
				FiberList[f].pointList.push_back(*p);

		}


	}

	void Resample(float spacing)
	{
		point3D<float> p0, p1;
		vector<point3D<float> > newPointList;
		for(unsigned int f=0; f<FiberList.size(); f++)
		{
			newPointList.clear();
			p0 = NodeList[FiberList[f].n0].p;
			for(unsigned int p=0; p<FiberList[f].pointList.size(); p++)
			{
				p1 = FiberList[f].pointList[p];
				if( (p1 - p0).Length() >= spacing )
				{
					newPointList.push_back(p1);
				}
			}
			FiberList[f].pointList = newPointList;
		}
	}
	//network comparison
	CoreGraphList CompareNetworks(rtsFiberNetwork* testNetwork, float sigmaG, float sigmaC, float &gFPR, float &gFNR, float &cFPR, float &cFNR)
	{
		//create point clouds that densely sample each network
		vector<geometryPoint> netPointList0, netPointList1;
		netPointList0 = getNetPointSamples(sigmaG);
		netPointList1 = testNetwork->getNetPointSamples(sigmaG);

		//compute the L1 distance between vertices in one network to the point cloud representing the other network
		KD_ComputeEnvelopeDistance(testNetwork, &netPointList0, sigmaG);
		KD_ComputeEnvelopeDistance(this, &netPointList1, sigmaG);

		//compute the geometry metric using the distance values for each vertex
		//float FPR, FNR;
		gFNR = GeometryMetric(this, sigmaG);
		gFPR = GeometryMetric(testNetwork, sigmaG);

		CoreGraphList core;
		core = NEW_ComputeTopology(testNetwork, sigmaC);

		//Changes made by James Burck (Thanks James!)--------
		float TP = (float)core.size();

		float TPandFP = (float)FiberList.size();      // formerly P, actaully TPandFN
		float TPandFN = (float)testNetwork->FiberList.size();   // actually TPandFP

		cFNR = (TPandFN - TP) / TPandFN;
		cFPR = (TPandFP - TP) / TPandFP;
		//---------------------------------------------------

		return core;
	}


};

void rtsFiberNetwork::initTopologyGraph(vector<topologyNode>* Nodes, vector<topologyEdge>* Edges, rtsFiberNetwork* network)
{
	/*This function constructs a graph based on the given network.*/
	Nodes->clear();
	Edges->clear();

	topologyNode node;
	topologyEdge edge;
	//for each node in the fiber network, construct a topologyNode
	for(unsigned int n=0; n<network->NodeList.size(); n++)
	{
		node.compatible = 0;
		node.label = RTS_TOPOLOGY_NODE_INVALID;
		node.p = network->NodeList[n].p;
		Nodes->push_back(node);
	}

	//now fill in all the edges
	for(unsigned int f=0; f<network->FiberList.size(); f++)
	{
		edge.n0 = network->FiberList[f].n0;
		edge.error = network->FiberList[f].error;
		edge.n1 = network->FiberList[f].n1;
		edge.label = RTS_TOPOLOGY_EDGE_EXIST;

		//attach the edge to each connected node in the node list
		(*Nodes)[edge.n0].connections.push_back(Edges->size());
		(*Nodes)[edge.n1].connections.push_back(Edges->size());

		//insert the edge into the list
		Edges->push_back(edge);
	}

}

void rtsFiberNetwork::BF_ComputeL1Distance(vector<geometryPoint>* N0, vector<geometryPoint>*N1)
{
	unsigned int i, j;
	vector3D<float> v;
	float dist;
	for(i=0; i<N0->size(); i++)
	{
		for(j=0; j<N1->size(); j++)
		{
			v = (*N0)[i].p - (*N1)[j].p;
			dist = v.Length();
			if(dist < (*N0)[i].dist)
				(*N0)[i].dist = dist;
			if(dist < (*N1)[j].dist)
				(*N1)[j].dist = dist;
		}
	}
}

void rtsFiberNetwork::BD_ComputeL1Distance(vector<geometryPoint>* N0, vector<geometryPoint>*N1)
{
	//build the point arrays
	ANNpointArray dataPts0 = annAllocPts(N0->size(), 3);
	for(unsigned int i=0; i<N0->size(); i++)
	{
		dataPts0[i][0] = (*N0)[i].p.x;
		dataPts0[i][1] = (*N0)[i].p.y;
		dataPts0[i][2] = (*N0)[i].p.z;
	}
	ANNpointArray dataPts1 = annAllocPts(N1->size(), 3);
	for(unsigned int i=0; i<N1->size(); i++)
	{
		dataPts1[i][0] = (*N1)[i].p.x;
		dataPts1[i][1] = (*N1)[i].p.y;
		dataPts1[i][2] = (*N1)[i].p.z;
	}

	//create ANN variables
	ANNbd_tree* bdTree;
	ANNpoint queryPt = annAllocPt(3);
	ANNidxArray nearestIdx = new ANNidx[1];
	ANNdistArray nearestDist = new ANNdist[1];

	//compare network 0 to network 1
	//bdTree = new ANNkd_tree(dataPts0, N0->size(), 3);
	bdTree = new ANNbd_tree(dataPts0, N0->size(), 3);
	for(unsigned int i=0; i<N1->size(); i++)
	{
		queryPt[0] = (*N1)[i].p.x;
		queryPt[1] = (*N1)[i].p.y;
		queryPt[2] = (*N1)[i].p.z;
		bdTree->annkSearch(queryPt, 1, nearestIdx, nearestDist);
		(*N1)[i].dist = sqrtf((float)nearestDist[0]);
	}
	delete bdTree;

	//compare network 1 to network 0
	bdTree = new ANNbd_tree(dataPts1, N1->size(), 3);
	for(unsigned int i=0; i<N1->size(); i++)
	{
		queryPt[0] = (*N0)[i].p.x;
		queryPt[1] = (*N0)[i].p.y;
		queryPt[2] = (*N0)[i].p.z;
		bdTree->annkSearch(queryPt, 1, nearestIdx, nearestDist);
		(*N0)[i].dist = sqrtf((float)nearestDist[0]);
	}
	delete bdTree;

	annClose();
}

void rtsFiberNetwork::MapDeviationToNetwork(vector<geometryPoint>* source)
{


}

void rtsFiberNetwork::topLabelNodes(vector<topologyNode>* N0, vector<topologyNode>* N1, float sigma)
{
	unsigned int i0, i1;
	vector3D<float> v;
	float min_d;
	unsigned int min_i;
	for(i0=0; i0 < N0->size(); i0++)
	{
		v = (*N0)[i0].p - (*N1)[0].p;
		min_d = v.Length();
		min_i = 0;
		for(i1=0; i1 < N1->size(); i1++)
		{
			v = (*N0)[i0].p - (*N1)[i1].p;
			if(v.Length() < min_d)
			{
				min_d = v.Length();
				min_i = i1;
			}
		}
		//if the minimum distance from point i0 is less than sigma
		if(min_d < sigma)
		{
			(*N0)[i0].label = RTS_TOPOLOGY_NODE_VALID;
			(*N0)[i0].compatible = min_i;
		}
	}
}

bool rtsFiberNetwork::topDetectEdge(vector<topologyNode>* NodeList, vector<topologyEdge>* EdgeList, unsigned int node0, unsigned int node1)
{
	//This function determines if there is an edge linking node0 and node1
	list<unsigned int>::iterator i;
	for(i = (*NodeList)[node0].connections.begin(); i!=(*NodeList)[node0].connections.end(); i++)
	{
		if( ((*EdgeList)[*i].n0 == node0 && (*EdgeList)[*i].n1 == node1) ||
			((*EdgeList)[*i].n0 == node1 && (*EdgeList)[*i].n1 == node0) )
			return true;
	}
	return false;
}
bool rtsFiberNetwork::topDeleteEdge(vector<topologyNode>* NodeList, vector<topologyEdge>* EdgeList, unsigned int node0, unsigned int node1)
{
	//this function deletes the first edge found linking node0 and node1
	list<unsigned int>::iterator i;
	unsigned int edge_id;
	for(i = (*NodeList)[node0].connections.begin(); i!=(*NodeList)[node0].connections.end(); i++)
	{
		if( (*EdgeList)[*i].n0 == node1 || (*EdgeList)[*i].n1 == node1 )
		{
			//delete the edge

			edge_id = *i;
			//remove the edge from node0 and node1
			(*NodeList)[node0].connections.remove(edge_id);
			(*NodeList)[node1].connections.remove(edge_id);
			//remove the edge
			(*EdgeList)[edge_id].label = RTS_TOPOLOGY_EDGE_NOEXIST;
			return true;
		}
	}
	return false;

}
bool rtsFiberNetwork::topMergeNode(vector<topologyNode>* NodeList, vector<topologyEdge>* EdgeList, unsigned int node)
{
	/*this function merges a specific node with it's neighbor based on the following rules:
	1) If the node is invalid, remove adjacent edge with the highest error
	2) If the node is valid, merge with adjacent compatible node, removing the edge with the largest error
	*/

	//if the node doesn't exist, just return.
	if( (*NodeList)[node].label == RTS_TOPOLOGY_NODE_NOEXIST) return false;
	//if this node isn't connected to anything, just remove the node
	if( (*NodeList)[node].connections.size() == 0)
	{
		(*NodeList)[node].label = RTS_TOPOLOGY_NODE_NOEXIST;
		return false;
	}

	//FIND THE DESTINATION NODE
	//create the destination node
	unsigned int edge_to_remove;

	//if the node is invalid, find the edge with the highest error
	if( (*NodeList)[node].label == RTS_TOPOLOGY_NODE_INVALID)
	{
		list<unsigned int>::iterator i;
		float highest_error = 0.0;
		for(i=(*NodeList)[node].connections.begin(); i!=(*NodeList)[node].connections.end(); i++)
		{
			//if the current edge has a higher error, record it
			if((*EdgeList)[(*i)].error >= highest_error)
			{
				highest_error = (*EdgeList)[(*i)].error;
				edge_to_remove = (*i);
			}
		}
	}

	//if the node is valid, find the compatible edge with the highest error
	if( (*NodeList)[node].label == RTS_TOPOLOGY_NODE_VALID)
	{
		list<unsigned int>::iterator i;
		float highest_error = 0.0;
		bool compatible_detected = false;
		unsigned int node0, node1;
		for(i=(*NodeList)[node].connections.begin(); i!=(*NodeList)[node].connections.end(); i++)
		{
			node0 = (*EdgeList)[(*i)].n0;
			node1 = (*EdgeList)[(*i)].n1;
			//find a compatible edge with the highest weight
			if((*NodeList)[node0].label == RTS_TOPOLOGY_NODE_VALID && (*NodeList)[node1].label == RTS_TOPOLOGY_NODE_VALID &&
				(*NodeList)[node0].compatible == (*NodeList)[node1].compatible && (*EdgeList)[(*i)].error >= highest_error)
			{
				highest_error = (*EdgeList)[(*i)].error;
				edge_to_remove = (*i);
				compatible_detected = true;
			}
		}
		//if a compatible node was not attached, just leave the node and return
		if(!compatible_detected) return false;
	}

	//PERFORM THE MERGE

	//find the node that we are merging to
	unsigned int merge_to;
	if((*EdgeList)[edge_to_remove].n0 == node)
		merge_to = (*EdgeList)[edge_to_remove].n1;
	else
		merge_to = (*EdgeList)[edge_to_remove].n0;

	list<unsigned int>::iterator i;
	//remove the edge from 'node'
	for(i = (*NodeList)[node].connections.begin(); i!=(*NodeList)[node].connections.end(); i++)
		if((*i) == edge_to_remove)
		{
			(*NodeList)[node].connections.erase(i);
			break;
		}
	//remove the edge from 'merge_to'
	for(i = (*NodeList)[merge_to].connections.begin(); i!=(*NodeList)[merge_to].connections.end(); i++)
		if((*i) == edge_to_remove)
		{
			(*NodeList)[merge_to].connections.erase(i);
			break;
		}

	//update all of the edges connected to 'node'
	for(i = (*NodeList)[node].connections.begin(); i!=(*NodeList)[node].connections.end(); i++)
	{
		if((*EdgeList)[(*i)].n0 == node)
			(*EdgeList)[(*i)].n0 = merge_to;
		else
			(*EdgeList)[(*i)].n1 = merge_to;
	}
	//add all edges in 'node' to the edges in 'merge_to'
	for(i = (*NodeList)[node].connections.begin(); i!=(*NodeList)[node].connections.end(); i++)
	{
		(*NodeList)[merge_to].connections.push_back( (*i) );
	}
	//sort the list and remove duplicates
	//duplicates occur if two merged points were connected by multiple edges
	(*NodeList)[merge_to].connections.sort();
	(*NodeList)[merge_to].connections.unique();

	//get rid of 'node'
	(*NodeList)[node].connections.clear();
	(*NodeList)[node].label = RTS_TOPOLOGY_NODE_NOEXIST;

	//remove the edge
	(*EdgeList)[edge_to_remove].label = RTS_TOPOLOGY_EDGE_NOEXIST;
	return true;
}
int rtsFiberNetwork::topCollapse(vector<topologyNode>* NodeList, vector<topologyEdge>*EdgeList)
{
	unsigned int n;
	unsigned int topology_changes = 0;
	unsigned int num_connections;
	bool node_merged = false;
	for(n=0; n<NodeList->size(); n++)
	{
		//if this node is the end of a barb, mark it as a topology change
		num_connections = (*NodeList)[n].connections.size();
		node_merged = topMergeNode(NodeList, EdgeList, n);
		if(num_connections == 1 && node_merged == true)
			topology_changes++;

	}

	return topology_changes;
}
void rtsFiberNetwork::MY_ComputeTopology(rtsFiberNetwork* testNetwork, float sigma)
{
	//initialize the topology graphs
	vector<topologyNode> GT_nodes;
	vector<topologyEdge> GT_edges;
	initTopologyGraph(&GT_nodes, &GT_edges, this);
	vector<topologyNode> T_nodes;
	vector<topologyEdge> T_edges;
	initTopologyGraph(&T_nodes, &T_edges, testNetwork);

	//label the nodes in each list as VALID or INVALID
	//this function also determines node compatibility in the Test array
	topLabelNodes(&GT_nodes, &T_nodes, sigma);
	topLabelNodes(&T_nodes, &GT_nodes, sigma);

	//copy the error to the fiber networks
	for(unsigned int i=0; i<T_nodes.size(); i++)
	{
		if(T_nodes[i].label == RTS_TOPOLOGY_NODE_VALID)
			testNetwork->NodeList[i].error = 0.0;
		else
			testNetwork->NodeList[i].error = 1.0;
	}
	for(unsigned int i=0; i<GT_nodes.size(); i++)
	{
		if(GT_nodes[i].label == RTS_TOPOLOGY_NODE_VALID)
			this->NodeList[i].error = 0.0;
		else
			this->NodeList[i].error = 1.0;
	}

	unsigned int FP_edges = topCollapse(&T_nodes, &T_edges);
	unsigned int FN_edges = topCollapse(&GT_nodes, &GT_edges);

	//mark all the nodes in T as invalid again
	//this will be used to compute topological errors
	unsigned int n;
	for(n=0; n<T_nodes.size(); n++)
	{
		T_nodes[n].label = RTS_TOPOLOGY_NODE_INVALID;
	}

	//for each edge in T
	unsigned int e;
	for(e=0; e<T_edges.size(); e++)
	{
		//if the edge exists in T
		if(T_edges[e].label == RTS_TOPOLOGY_EDGE_EXIST)
		{
			//if the the edge exists in GT
			if(topDetectEdge(&GT_nodes, &GT_edges, T_nodes[T_edges[e].n0].compatible, T_nodes[T_edges[e].n1].compatible))
			{
				//set both nodes to valid, in order to detect topological errors
				T_nodes[T_edges[e].n0].label = RTS_TOPOLOGY_NODE_VALID;
				T_nodes[T_edges[e].n1].label = RTS_TOPOLOGY_NODE_VALID;
				//delete the edge in GT
				topDeleteEdge(&GT_nodes, &GT_edges, T_nodes[T_edges[e].n0].compatible, T_nodes[T_edges[e].n1].compatible);
			}
			else
				FP_edges++;
		}
	}

	//run through all edges in GT
	//for each edge that still exists, increment the FN_edges counter
	for(e=0; e<GT_edges.size(); e++)
	{
		if(GT_edges[e].label == RTS_TOPOLOGY_EDGE_EXIST)
			FN_edges++;
	}

	//find topological errors
	//these are nodes that are valid AND are duplicates
	list<unsigned int> validPoints;
	for(n=0; n<T_nodes.size(); n++)
		if(T_nodes[n].label == RTS_TOPOLOGY_NODE_VALID)
			validPoints.push_back(T_nodes[n].compatible);
	validPoints.sort();
	list<unsigned int>::iterator i;
	unsigned int last;
	unsigned int topErrors = 0;
	for(i = validPoints.begin(); i != validPoints.end(); i++)
	{
		if(i != validPoints.begin())
			if(*i == last)
				topErrors++;
		last = *i;
	}



	cout<<"False Positive Edges: "<<FP_edges<<endl;
	cout<<"False Negative Edges: "<<FN_edges<<endl;
	cout<<"Topological Errors: "<<topErrors<<endl;


}
TopologyGraph rtsFiberNetwork::BOOST_RemapGraph(TopologyGraph& G, rtsFiberNetwork* network)
{
	//create the graph
	TopologyGraph result(network->NodeList.size());

	//add all edges in G to result (based on compatibility (original index), NOT current index)
	graph_traits<TopologyGraph>::edge_iterator ei, ei_end;
	graph_traits<TopologyGraph>::vertex_descriptor v0, v1;
	int idx0, idx1;
	for(tie(ei, ei_end) = edges(G); ei != ei_end; ei++)
	{
		v0 = source(*ei, G);
		v1 = target(*ei, G);
		idx0 = get(vertex_color_t(), G, v0);
		idx1 = get(vertex_color_t(), G, v1);
		add_edge(idx0, idx1, result);
	}
	return result;

}
point3D<float> rtsFiberNetwork::getNodeCoord(int fiber, bool node)
{
	//Return the node coordinate based on a fiber (edge)
	//node value is false = source, true = dest

	if(node)
		return getNodeCoord(FiberList[fiber].n1);
	else
		return getNodeCoord(FiberList[fiber].n0);

}

double rtsFiberNetwork::getTotalLength()
{
	//go through each fiber and measure the length
	int num_fibers = FiberList.size();
	double length = 0.0;
	int f;
	for(f=0; f<num_fibers; f++)
		length += getFiberLength(f);
	return length;
}
point3D<float> rtsFiberNetwork::getFiberPoint(unsigned int fiber, unsigned int point)
{
	if(point == 0)
		return NodeList[FiberList[fiber].n0].p;
	if(point == FiberList[fiber].pointList.size() + 1)
		return NodeList[FiberList[fiber].n1].p;

	return FiberList[fiber].pointList[point-1];

}

void rtsFiberNetwork::ComputeBoundingVolume()
{
	//find the bounding volume for the nodes
	min_pos = NodeList[0].p;
	max_pos = NodeList[0].p;
	for(unsigned int n=0; n<NodeList.size(); n++)
	{
		if(NodeList[n].p.x < min_pos.x)
			min_pos.x = NodeList[n].p.x;
		if(NodeList[n].p.y < min_pos.y)
			min_pos.y = NodeList[n].p.y;
		if(NodeList[n].p.z < min_pos.z)
			min_pos.z = NodeList[n].p.z;

		if(NodeList[n].p.x > max_pos.x)
			max_pos.x = NodeList[n].p.x;
		if(NodeList[n].p.y > max_pos.y)
			max_pos.y = NodeList[n].p.y;
		if(NodeList[n].p.z > max_pos.z)
			max_pos.z = NodeList[n].p.z;
	}

	//combine with the bounding volume for the fibers
	for(unsigned int f=0; f<FiberList.size(); f++)
	{
		for(unsigned int p=0; p<FiberList[f].pointList.size(); p++)
		{
			if(FiberList[f].pointList[p].x < min_pos.x)
				min_pos.x = FiberList[f].pointList[p].x;
			if(FiberList[f].pointList[p].y < min_pos.y)
				min_pos.y = FiberList[f].pointList[p].y;
			if(FiberList[f].pointList[p].z < min_pos.z)
				min_pos.z = FiberList[f].pointList[p].z;

			if(FiberList[f].pointList[p].x > max_pos.x)
				max_pos.x = FiberList[f].pointList[p].x;
			if(FiberList[f].pointList[p].y > max_pos.y)
				max_pos.y = FiberList[f].pointList[p].y;
			if(FiberList[f].pointList[p].z > max_pos.z)
				max_pos.z = FiberList[f].pointList[p].z;
		}
	}

}
void rtsFiberNetwork::Translate(float x, float y, float z)
{
	//translates the network while maintaining all connectivity

	//create a translation vector
	vector3D<float> translate(x, y, z);

	//translate all nodes
	int num_nodes = NodeList.size();
	int n;
	for(n=0; n<num_nodes; n++)
		NodeList[n].p = NodeList[n].p + translate;

	int num_fibers = FiberList.size();
	int num_points;
	int f, p;
	for(f=0; f<num_fibers; f++)
	{
		num_points = FiberList[f].pointList.size();
		for(p=0; p<num_points; p++)
			FiberList[f].pointList[p] = FiberList[f].pointList[p] + translate;
	}

	//translate the bounding box
	min_pos = min_pos + translate;
	max_pos = max_pos + translate;

}

#endif
