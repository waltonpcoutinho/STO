#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <vector>
#include <time.h>
#include <sys/time.h>
#include <climits>
#include <cfloat>
#include <random>
#include <functional>

#include "Data.h"
#include "Solution.h"

#define INF INT_MAX
#define INFF DBL_MAX
#define epsl 0.001
#define penal 1
#define maxtb 10000

#define DEBUG 0

//macro for pausing, showing file name and line number
#define PAUSE do {                                      \
   printf("paused at %s, line %d\n", __FILE__, __LINE__); \
   getchar();                                           \
}while(0)

using namespace std;

//forward declaration of class SP.h
//because SP includes Util and Util 
//includes SP (compiler gets lost)
//class SP;

//function declarations
void setSeed();
bool randomBool();
double fRand(double, double);
void randomSeq(int*, int);
double cpuTime();
double wallTime();
double l2_norm(const double, const double, const double);
template<typename Type> void printMatrix(Type, std::size_t, std::size_t);
void printVec(double const*,int, string);
void printVec(int const*,int, string);
void printVec(vector<int> vec, string text);
void printRoutes(vector<glider> gliders);
double routeTotalFlightTime(glider uav);
double routeTotalError(glider uav);
double routeAvgStep(glider uav);
void writeRes2File(Data const*, Solution const*);
void writeInfList(Data const*);

void writeSol2File(Data const*, vector<glider>);


#endif
