/*
-- This file contains single-source and recursive single pair rolesim* similarity functions.This version:
---- uses a Threshold (theta) to stop computing similarity after fixed number of iterations (itt>2)
	 for pairs whose (similarity in current iteration)-(similarity in previous iteration)<theta
---- uses a list of equivalent nodes (equiSet) to ignore duplicate similarity computations for nodes
	 that their similarity is always equal.
---- uses a threshold (delta) to control the memory complexity by storing the similarity of pairs whose
	 (outdegree>delta)
---- in order to record the similarities this functions use a nested-dictionary (unordered_map smap)
	 that for each pair (a,b) records a_id,b_id,tuple(ittr,haspruned,sim)
	 - recording ittr is the itteration that the similarity is computed in order to ignore computing similarities
	   in lower itterations
	 - recording haspruned(T/F) is to ignore computing similarity if it is pruned by the theta threshold
*/

#include "stdafx.h"

void SSRSStar_OptV2_Trsh(PNGraph& G, vector<int>& d, int k, float beta, float lambda, int a, int b, double& s, unordered_map<int, unordered_map<int, tuple<int, bool, double>>>& smap, long& count, double theta, int delta, bool lastLevel, vector<int>& equiSet, vector<long>& prunedCount);
void SSRSStar_OptV2_Trsh(PNGraph& G, vector<int>& degVec, const int kmax, float beta, float lambda, int q, vector<double>& sq, long& count, double theta, int delta, bool lastLevel, vector<int>& equiSet, vector<long>& prunedCount)
{
	auto n = degVec.size();
	double s(0.0);
	clock_t t1, t2;
	pair<int, double> tempPair;
	unordered_map<int, unordered_map<int, tuple<int, bool, double>>> smap;
	for (TNGraph::TNodeI a = G->BegNI(); a < G->EndNI(); a++)
	{
		t1 = clock();
		const int a_id = a.GetId();
		if (a_id < n * 0.4)
		{
			if (degVec[a_id] == 0 || degVec[q] == 0 || kmax <= 0) {
				s = (1.0 - beta);
			}
			else
			{
				int a_id_new = equiSet[a_id];
				int q_new = equiSet[q];
				if (a_id_new > q_new)
					swap(a_id_new, q_new);
				SSRSStar_OptV2_Trsh(G, degVec, kmax, beta, lambda, a_id_new, q_new, s, smap, count, theta, delta, lastLevel, equiSet, prunedCount);
			}
			//	sq.push_back(s);
			sq[a_id] = s;
		}
		t2 = clock();
		//cout << a_id << "\t" << s << "\t" << (t2 - t1) / (float)CLOCKS_PER_SEC << endl;
	}
	cout << q << " ";
}
void SSRSStar_OptV2_Trsh(PNGraph& G, vector<int>& degVec, int k, float beta, float lambda, int a, int b, double& s, unordered_map<int, unordered_map<int, tuple<int, bool, double>>>& smap, long& count, double theta, int delta, bool lastLevel, vector<int>& equiSet, vector<long>& prunedCount)
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
				a_nbr = equiSet[a_nbr];
				b_nbr = equiSet[b_nbr];
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
					SSRSStar_OptV2_Trsh(G, degVec, k - 1, beta, lambda, n1, n2, sim, smap, count, theta, delta, lastLevel, equiSet, prunedCount);
				}
				else
				{
					auto it1 = it->second.find(n2);
					if (it1 == it->second.end())
						SSRSStar_OptV2_Trsh(G, degVec, k - 1, beta, lambda, n1, n2, sim, smap, count, theta, delta, lastLevel, equiSet, prunedCount);
					else
					{
						int  prek = std::get<0>(it1->second);
						bool flag = std::get<1>(it1->second);
						if (prek < k - 1 && !flag)
						{
							SSRSStar_OptV2_Trsh(G, degVec, k - 1, beta, lambda, n1, n2, sim, smap, count, theta, delta, lastLevel, equiSet, prunedCount);
						}
						else
						{
							sim = std::get<2>(it1->second);
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
	//cout <<"curr:"<< s;
	bool f = false;
	tuple<int, bool, double> value;
	auto it = smap.find(a);
	if (it != smap.end())
	{
		auto it1 = it->second.find(b);
		if (it1 != it->second.end())
		{
			double preSim = get<2>(it1->second);
			//cout << "\t pre:" << preSim;
			double diff = abs(preSim - s);//(10a)
			/*if (diff > 0)
				cout << diff << "\t" << preSim << "\t" <<s << endl;*/
			if (diff < theta && k >= 2)
			{
				f = true;
			}
			if (s < ((1.0 - beta) + theta)) //(10b)
			{
				f = true;//flag
				s = 1.0 - beta;//sim
			}
		}
	}
	if (f)
		prunedCount[k]++;
	//cout << endl;
	value = make_tuple(k, f, s);
	if (aa.GetOutDeg() * bb.GetOutDeg() > delta)
	{
		smap[a][b] = value;
	}
	count++;
}

