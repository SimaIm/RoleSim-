/*
-- This file contains threshold-based rolesim* similarity function.This version:
   --uses a Threshold (theta) to stop computing similarity after fixed number of iterations (itt>2)
	 for pairs whose (similarity in current iteration)-(similarity in previous iteration)<theta

*/
#include "stdafx.h"


void RoleSimStar_Trsh(PNGraph& G, vector<vector<double>>& SimMtx, int kmax, float beta, float lambda,float theta)
{
	vector<int> assignment;
	int nodes_count = G->GetNodes();
	vector<int> pruned_count(kmax);
	vector<vector<double>> PreSimMtx;
	vector< vector<bool> > ThrshMtx(nodes_count, vector<bool >(nodes_count,true));
	for (int k = 0; k < kmax; k++)
	{
		PreSimMtx = SimMtx;
		for (TNGraph::TNodeI a = G->BegNI(); a < G->EndNI(); a++)
		{
			const int a_id = a.GetId();
			const int a_ideg = a.GetInDeg();
			for (TNGraph::TNodeI b = a; b < G->EndNI(); b++)
			{
				double sum = 0.0, partA = 0.0, partB = 0.0;
				double maxMatch = 0;
				const int b_id = b.GetId();
				const int b_ideg = b.GetInDeg();
				if (ThrshMtx[a_id][b_id])
				{
					if (a_ideg == 0 || b_ideg == 0)
					{
						SimMtx[a_id][b_id] = (1.0 - beta);
						SimMtx[b_id][a_id] = (1.0 - beta);
					}
					else
					{
						if (a_ideg == 1 && b_ideg == 1)
						{
							int a_nbr = a.GetInNId(0);
							int b_nbr = b.GetInNId(0);
							maxMatch = PreSimMtx[a_id][b_id];
						}
						else
						{
							TNGraph::TNodeI tn1, tn2;
							int64_t tn1_ideg, tn2_ideg;
							if (a_ideg <= b_ideg)
							{
								tn1 = a;
								tn1_ideg = a_ideg;
								tn2 = b;
								tn2_ideg = b_ideg;
							}
							else
							{
								tn1 = b;
								tn1_ideg = b_ideg;
								tn2 = a;
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
									nbrMat[bn + (tn2_ideg * an)] = PreSimMtx[a_nbr][b_nbr] * (-1);
									sum = sum + PreSimMtx[a_nbr][b_nbr];
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
						SimMtx[a_id][b_id] = beta * (partA + partB) + (1.0 - beta);
						SimMtx[b_id][a_id] = SimMtx[a_id][b_id];
					}
					if (k >= 2 && (PreSimMtx[a_id][b_id] - SimMtx[a_id][b_id]) <= theta)
					{
						pruned_count[k]++;
						ThrshMtx[a_id][b_id] = false;
						//ThrshMtx[b_id][a_id] = false;
						if (SimMtx[a_id][b_id] <= ((1 - beta) + theta))
						{
							SimMtx[a_id][b_id] = 1 - beta;
							SimMtx[b_id][a_id] = 1 - beta;
						}
					}
					else if (SimMtx[a_id][b_id] <= ((1 - beta) + theta))
					{
						pruned_count[k]++;
						ThrshMtx[a_id][b_id] = false;
						SimMtx[a_id][b_id] = 1 - beta;
						SimMtx[b_id][a_id] = 1 - beta;
					}
				}
			}// b
		}//a
	}//k
}