#pragma once
#ifndef FUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FUNCTIONS_H
#include "stdafx.h"

//================== ALL-PAIR SIMILARITY FUNCTIONS 
void RoleSim(PNGraph& G, vector<vector<double>>& SimMtx, int kmax, float beta);
void RoleSimStar(PNGraph& G, vector<vector<double>>& SimMtx, int kmax, float beta, float lambda);
void RoleSimStar_Trsh(PNGraph& G, vector<vector<double>>& SimMtx, int kmax, float beta, float lambda, float theta);


//================== SINGLE-SOURCE SIMILARITY FUNCTIONS 
void SSRSStar_OptV1(PNGraph& G, vector<int>& degVec, int kmax, float beta, float lambda, int q, vector<double>& sq, long& count, int theta, bool lastLevel);
void SSRSStar_OptV2(PNGraph& G, vector<int>& degVec, int kmax, float beta, float lambda, int q, vector<double>& sq, long& count, int theta, bool lastLevel, vector<int>& equiSet);
void SSRSStar_OptV2_Trsh(PNGraph& G, vector<int>& degVec, const int kmax, float beta, float lambda, int q, vector<double>& sq, long& count, double theta, int delta, bool lastLevel, vector<int>& equiSet, vector<long>& prunedCount);



//================== GENERAL FUNCTIONS 
void BuildGraph      (PNGraph& G, vector<int>& indegVec);
void FindEqueivalents(PNGraph& G, vector<int>& equiList);
#endif