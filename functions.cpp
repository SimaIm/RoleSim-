#include "stdafx.h"
#include <fstream>
#include "CTrie.h"
// A Trie node

void BuildGraph(PNGraph& G, vector<int>& indegVec)
{
	int count = 0;
	for (TNGraph::TNodeI v = G->BegNI(); v < G->EndNI(); v++) {
		const int v_id = v.GetId();
		clock_t t1, t2, t3;
		const int v_deg = v.GetInDeg();
		indegVec[v_id] = v_deg;
	}
	//cout << count;
}

void FindEqueivalents(PNGraph& G, vector<int>& equiSet)
{
	//FILE* Fp;
	//string Path = "E:/PhD_Research/Implementations/Data/RealData/email-dnc_Equi.txt";
	//Fp = fopen(Path.c_str(), "w");
	equiSet.resize(G->GetNodes());
	pair<int, vector<int>> max;
	vector<int> temp(0);
	max = make_pair(0, temp);
	int count = 0;
	CTrie* head = new CTrie();
	for (TNGraph::TNodeI v = G->BegNI(); v < G->EndNI(); v++)
	{
		const int v_id = v.GetId();
		const int v_deg = v.GetInDeg();
		if (v_deg)
		{
			vector<int>  nbrSet(v_deg);
			for (int e = 0; e < v_deg; e++) {
				nbrSet[e] = v.GetInNId(e);
			}
			max = (v_deg > max.first) ? make_pair(v_deg, nbrSet) : max;
			int eqID = -1;
			bool exist = head->search(nbrSet, eqID);
			if (!exist)
			{
				equiSet[v_id] = v_id;
				head->insert(nbrSet, v_id);
			}
			else
			{
				equiSet[v_id] = eqID;
				count++;
			}
		}
		else
		{
			equiSet[v_id] = v_id;
		}
	}
	delete head;
	//fclose(Fp);
	//cout << "-Number of duplicates " << count << endl;
}
