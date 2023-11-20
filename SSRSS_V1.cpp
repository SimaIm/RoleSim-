/*
-- This file contains single-source and recursive single pair rolesim* similarity functions.This version:
---- uses a threshold (delta) to control the memory complexity by storing the similarity of pairs whose
	 (outdegree>delta)
---- in order to record the similarities this functions use a nested-dictionary (unordered_map smap)
	 that for each pair (a,b) records a_id,b_id,tuple(ittr,haspruned,sim)
	 - recording ittr is the itteration that the similarity is computed in order to ignore computing similarities
	   in lower itterations
*/

#include "stdafx.h"
void SPRSStar_OptV1_Rec(PNGraph& G, vector<int>& d, int k, float beta, float lambda, int a, int b,
	double& s, unordered_map<int, unordered_map<int, pair<int, double>>>& smap,
	long& count, int delta, bool lastLevel);
void SSRSStar_OptV1(PNGraph& G, vector<int>& degVec, int kmax, float beta, float lambda, int q,
	vector<double>& sq, long& count, int delta, bool lastLevel)
{
	auto n = degVec.size();
	double s(0.0);
	clock_t t1, t2;
	pair<int, double> tempPair;
	unordered_map<int, unordered_map<int, pair<int, double>>> smap;
	for (TNGraph::TNodeI a = G->BegNI(); a < G->EndNI(); a++)
	{
		t1 = clock();
		const int a_id = a.GetId();
		if (degVec[a_id] == 0 || degVec[q] == 0 || kmax <= 0) {
			s = (1.0 - beta);
		}
		else {
			int a_id_new = a_id;
			int q_new = q;
			if (a_id_new > q_new)
				swap(a_id_new, q_new);
			SPRSStar_OptV1_Rec(G, degVec, kmax, beta, lambda, a_id_new, q_new, s, smap, count, delta, lastLevel);
		}
		sq[a_id] = s;
		t2 = clock();
		//cout << a_id << "\t" << degVec[a_id] << "\t" << (t2 - t1) / (float)CLOCKS_PER_SEC << endl;
	}
	//cout << endl;
}
void SPRSStar_OptV1_Rec(PNGraph& G, vector<int>& degVec, int k, float beta, float lambda, int a, int b,
	double& s, unordered_map<int, unordered_map<int, pair<int, double>>>& smap,
	long& count, int delta, bool lastLevel)
{
	int a_ideg = degVec[a];
	int b_ideg = degVec[b];
	if (a_ideg == 0 || b_ideg == 0 || k <= 0) {
		s = (1.0 - beta);
		return;
	}
	TNGraph::TNodeI aa = G->GetNI(a);
	TNGraph::TNodeI bb = G->GetNI(b);
	double maxMatch = 0;
	TNGraph::TNodeI tn1, tn2;
	int64_t tn1_ideg, tn2_ideg;
	double sum = 0.0, partA = 0.0, partB = 0.0;
	if ((k - 1) == 0 && lastLevel)
	{
		maxMatch = min(b_ideg, a_ideg) * (1.0 - beta);
		sum = (b_ideg * a_ideg) * (1.0 - beta);
	}
	else
	{
		if (a_ideg <= b_ideg)
		{
			tn1 = aa;
			tn1_ideg = a_ideg;
			tn2 = bb;
			tn2_ideg = b_ideg;
		}
		else
		{
			tn1 = bb;
			tn1_ideg = b_ideg;
			tn2 = aa;
			tn2_ideg = a_ideg;
		}
		double* nbrMat;
		nbrMat = new double[tn1_ideg * tn2_ideg];
		int64_t* col4row = new int64_t[tn1_ideg];
		for (int an = 0; an < tn1_ideg; an++)
		{
			int a_nbr = tn1.GetInNId(an);
			for (int bn = 0; bn < tn2_ideg; bn++)
			{
				int b_nbr = tn2.GetInNId(bn);
				int n1 = 0, n2 = 0;
				if (a_nbr <= b_nbr)
				{
					n1 = a_nbr;
					n2 = b_nbr;
				}
				else
				{
					n1 = b_nbr;
					n2 = a_nbr;
				}
				double sim(0.0);
				auto it = smap.find(n1);
				if (it == smap.end())
				{
					SPRSStar_OptV1_Rec(G, degVec, k - 1, beta, lambda, n1, n2, sim, smap, count, delta, lastLevel);
				}
				else
				{
					auto it1 = it->second.find(n2);
					if (it1 == it->second.end())
						SPRSStar_OptV1_Rec(G, degVec, k - 1, beta, lambda, n1, n2, sim, smap, count, delta, lastLevel);
					else
					{
						int prek = it1->second.first;
						if (prek != k - 1)
							SPRSStar_OptV1_Rec(G, degVec, k - 1, beta, lambda, n1, n2, sim, smap, count, delta, lastLevel);
						else
						{
							sim = it1->second.second;
						}
					}
				}
				nbrMat[bn + (tn2_ideg * an)] = sim * (-1.0);
				sum = sum + sim;
			}
		}
		solve_rectangular_linear_sum_assignment(tn1_ideg, tn2_ideg, nbrMat, col4row);
		for (int64_t i = 0; i < tn1_ideg; i++)
		{
			maxMatch = maxMatch + (nbrMat[col4row[i] + (tn2_ideg * i)] * (-1.0));
		}
		delete[] nbrMat;
		delete[] col4row;
	}
	partA = lambda * (maxMatch / (double)max(a_ideg, b_ideg));
	if (a_ideg == 1 && b_ideg == 1)
		partB = 0;
	else
		partB = (1.0 - lambda) * ((sum - maxMatch) / double((a_ideg * b_ideg) - min(a_ideg, b_ideg)));

	s = beta * (partA + partB) + (1.0 - beta);

	pair<int, double> value;
	value = make_pair(k, s);
	if (aa.GetOutDeg() * bb.GetOutDeg() >= delta)
		smap[a][b] = value;
	count++;
	//cout << k <<"-> "<< a << "," << b << "="<< value.second<<endl;
}
