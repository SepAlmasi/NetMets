// add the following to a cpp file:
// PerformanceData PD;

#pragma once

enum PerformanceDataType
{
	PD_DISPLAY=0,
	PD_SPS,
	PD_UNUSED0,

	//my stuff
	SUBDIVIDE_NETWORKS,
	BRUTE_FORCE_DISTANCE,
	LOG_N_DIST_BUILD0,
	LOG_N_DIST_SEARCH0,
	LOG_N_DIST_BUILD1,
	LOG_N_DIST_SEARCH1,
	LOG_N_DISTANCE,
	MY_TOPOLOGY,
	BOOST_TOPOLOGY,
	

	//end my stuff
	PERFORMANCE_DATA_TYPE_COUNT
};

static char PDTypeNames[][255] = {
		"Display            ",
		"Simulation Total   ",
		" ----------------- ",
		//my stuff
		"Network Subdivision",
		"Brute Force Dist.--",
		"    Build 0        ",
		"    Search 0       ",
		"    Build 1        ",
		"    Search 1       ",
		"ANN Dist.----------",
		"MY Topology--------",
		"BOOST Topology-----",
		//end my stuff
		
};

//-------------------------------------------------------------------------------

//#ifdef DISPLAY
//	#define PERFORMANCE_DATA_MEASUE_ENABLE
//#endif

//#ifdef PERFORMANCE_DATA_MEASUE_ENABLE

//-------------------------------------------------------------------------------

#include <stdio.h>
#include <windows.h>
#include <float.h>

#include <iostream>
#include <iomanip>

//-------------------------------------------------------------------------------

class PerformanceData
{
public:
	PerformanceData() { ClearAll(); QueryPerformanceFrequency(&cps); }
	~PerformanceData(){}

	void ClearAll()
	{
		for ( int i=0; i<PERFORMANCE_DATA_TYPE_COUNT; i++ ) {
			for ( int j=0; j<256; j++ ) times[i][j] = 0;
			pos[i] = 0;
			minTime[i] = 0xFFFFFFFF;
			maxTime[i] = 0;
			totalTime[i] = 0;
			dataReady[i] = false;
		}
	}

	void StartTimer( int type ) { QueryPerformanceCounter( &startTime[type] ); /*startTime[type] = GetTickCount();*/ }
	void EndTimer( int type ) {
		LARGE_INTEGER endTime;
		QueryPerformanceCounter( &endTime );
		double t = endTime.QuadPart - startTime[type].QuadPart;
		//unsigned int t = GetTickCount() - startTime[type];
		if ( t < minTime[type] ) minTime[type] = t;
		if ( t > maxTime[type] ) maxTime[type] = t;
		totalTime[type] -= times[type][ pos[type] ];
		times[type][ pos[type] ] = t;
		totalTime[type] += t;
		pos[type]++;
		if ( pos[type] == 0 ) dataReady[type] = true;
	}

	void PrintResult( ostream &os,int i=PERFORMANCE_DATA_TYPE_COUNT)
	{
		os.setf(ios::fixed);
		if ((i<PERFORMANCE_DATA_TYPE_COUNT)&&(i>=0)){
			double a = GetAvrgTime(i);
			if ( a ) 
				os<< PDTypeNames[i]<<" :  avrg="<<setw(8)<<setprecision(3)<<a<<"\tmin="<<setw(8)<<setprecision(3)<< GetMinTime(i) <<"\tmax="<<setw(8)<<setprecision(3)<< GetMaxTime(i) <<endl ;
			else 
				os<< PDTypeNames[i]<<" :  avrg=   -----\tmin=   -----\tmax=   -----"<<endl;
		}
	}

	void PrintResults( ostream &os)
	{
		for ( int i=0; i<PERFORMANCE_DATA_TYPE_COUNT; i++ ) 
			PrintResult(os,i);
	}

	double	GetLastTime( int type ) { return times[type][pos[type]]; }
	double	GetAvrgTime( int type ) { double a = 1000.0 * totalTime[type] / (float)cps.QuadPart / ( (dataReady[type]) ? 256.0 : (double)pos[type] ); return (_finite(a))? a:0; }
	double	GetMinTime( int type ) { return 1000.0 * minTime[type] / (float)cps.LowPart; }
	double	GetMaxTime( int type ) { return 1000.0 * maxTime[type] / (float)cps.LowPart; }

private:
	double times[PERFORMANCE_DATA_TYPE_COUNT][256];
	unsigned char pos[PERFORMANCE_DATA_TYPE_COUNT];
	LARGE_INTEGER startTime[PERFORMANCE_DATA_TYPE_COUNT];
	double minTime[ PERFORMANCE_DATA_TYPE_COUNT ];
	double maxTime[ PERFORMANCE_DATA_TYPE_COUNT ];
	double totalTime[ PERFORMANCE_DATA_TYPE_COUNT ];
	bool dataReady[ PERFORMANCE_DATA_TYPE_COUNT ];
	LARGE_INTEGER cps;
};

//-------------------------------------------------------------------------------
/*#else	

class PerformanceData{
public:
	PerformanceData() {;};
	~PerformanceData(){;};
	void ClearAll(){;};
	void StartTimer( int type ) {;};
	void EndTimer( int type ) {;};
	void PrintResults( ostream &os){;};
	void PrintResult( ostream &os, int i=PERFORMANCE_DATA_TYPE_COUNT){;};
	double	GetLastTime( int type ) { return 0.0; };
	double	GetAvrgTime( int type ) { return 0.0; };
	double	GetMinTime( int type ) { return 0.0; };
	double	GetMaxTime( int type ) { return 0.0; };
};

#endif
//-------------------------------------------------------------------------------
*/
extern PerformanceData PD;

//-------------------------------------------------------------------------------
