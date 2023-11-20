#include "stdafx.h"

int main(int argc, char* argv[]) {

    /* ------ Parameters---------------*/
    int kmax = 10;
    double beta = 0.8;// Damping factor
    double lambda = 0.7;// RS* lambda
    double delta = 0;// Threshold factor
    double theta = 0.05;// Memory Threshold factor
    long count = 0; // number of computed pairs using single-source
    vector<long> pruned_count; // number of pruned pairs using threshold per k
    double sim0 = 0.2;

    /* ------ Input Data---------------*/
    TStr fName = "ca-CSphd";//Cit-HepPh
    clock_t t0, t1;
    TStr fullname("E:/PhD_Research/Implementations/Data/RealData/" + fName + ".txt");

    cout << "> Loading Graph to SNAP ... " << endl;
    TStrHash< TInt > NodeNames;
    PNGraph G = TSnap::LoadEdgeListStr<PNGraph>(fullname, 0, 1, NodeNames);

    /* ------ Graph Info---------------*/
    int nodes = G->GetNodes();
    int edges = G->GetEdges();
    cout << "\tGraph Name = \t" << fName.CStr() << "\n"
        << "\t# of Nodes = \t" << nodes << "\n"
        << "\t# of Edges = \t" << edges << "\n"
        << "\tAve Degree = \t" << (float)edges / nodes << endl;

    /* ------ Build Graph*/
    cout << " Building Graph ..." << endl;
    vector<int> indegVec(nodes);
    BuildGraph(G, indegVec);
   
   /* ------ Finding Equivalents Nodes, IF USING VERSION2 */
    cout << " Finding Equivalents Nodes ..." << endl;
    vector<int> equiSet;
    FindEqueivalents(G, equiSet);

    /* ------ All-pairs Similarity computation ---------*/
    vector<vector<double>> simMtx(nodes, vector<double>(nodes, 0.2));//matrix similarity vector
    RoleSimStar(G, simMtx, kmax, beta, lambda);
    RoleSimStar_Trsh(G, simMtx, kmax, beta, lambda,theta);
    /* ------ Single_Source Similarity computation ---------*/
    int query = 1;
    vector<double> qSimVec(nodes, sim0);//single source similarity vector
    SSRSStar_OptV1( G, indegVec, kmax, beta, lambda, query, qSimVec, count, theta, true);
    SSRSStar_OptV2( G, indegVec, kmax, beta, lambda, query, qSimVec, count, theta, true, equiSet);
    SSRSStar_OptV2_Trsh(G, indegVec, kmax, beta, lambda, query, qSimVec, count, theta, delta, true, equiSet, pruned_count);
    
    return 0;
}