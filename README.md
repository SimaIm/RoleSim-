# RoleSim*

RoleSim* is a similarity measure introduced in the paper ``RoleSim*: Scaling axiomatic role-based similarity ranking on large graphs''. This measure combines the strengths of both SimRank and RoleSim while effectively addressing some of their limitations.
This repository contains the C++ implementation of the following algorithms:

- Variants of RoleSim* similarity measure, including the threshold-based and single-source versions.
- Original [SimRank](https://dl.acm.org/doi/10.1145/775047.775126) similarity measure.
- Original [RoleSim](https://dl.acm.org/doi/abs/10.1145/2020408.2020561) similarity measure.

## Paper Reference
For a detailed understanding of the methodology behind variants of RoleSim*, refer to the [Original Paper](https://link.springer.com/article/10.1007/s11280-021-00925-z).

## Configuration
All algorithms are implemented in C++. Variants of the RoleSim* algorithm use the Stanford Network Analysis Platform (SNAP), which provides general-purpose network analysis functions.
To use the RoleSim* variants with SNAP, it is necessary to configure Visual Studio. Follow the instructions in [this link](https://snap.stanford.edu/) for more details.


## List of Functions
RoleSimStar: all-pairs RS* similarity function.
SSRSStar_OptV1_Rec: Single-Source RS* similarity function
SSRSStar_OptV2_Rec: Single-Source RS* similarity function using the equivalent nodes to ignore repeated computations for pairs with equal neighbours.
SSRSStar_OptV2_Trsh: Threshold-based Single-Source RS* similarity function using the equivalent nodes to ignore repeated computations for pairs with equal neighbours and considering a threshold to stop computing similarity after a fixed number of iterations.    
FindEqueivalents: This function finds nodes with the equal in-neighbours.


## Usage
There are some initializations before calling given similarity functions.
First, we initialize parameters of the algorithm named according to the paper, and then loading the graph file using SNAP Library function. Third, computing an array of indegree’s for the input graphs.

  /* ------ Initializing Parameters---------------*/
    int kmax = 10; // # of iterations in the all-pairs version or depth in single-source version
    double beta = 0.8;// Damping factor
    double lambda = 0.7;// RS* lambda
    double delta = 0;// Threshold factor
    double theta = 0.05;// Memory Threshold factor

    /* ------ Load Network Data---------------*/
    TStrHash< TInt > NodeNames;
    PNGraph G = TSnap::LoadEdgeListStr<PNGraph>("NetworkFile.txt", 0, 1, NodeNames);

    /* ------ Extract Indegree Vector*/
    vector<int> indegVec(N);
    BuildGraph(G, indegVec);
   
To compute similarity of all pairs of nodes in the graph, we use RoleSimStar function from RSS.cpp file. The output of this function is a NxN matrix of similarities, where N is the number of nodes in the graph.

    /* ------ All-pairs Similarity computation ---------*/
    vector<vector<double>> simMtx(N, vector<double>(N, 0.2));// similarity matrix 
    RoleSimStar(G, simMtx, kmax, beta, lambda);


To compute similarity of a query node with all other nodes in the graph, we use SSRSStar_OptV1 (SSRSS_V1.cpp), SSRSStar_OptV2 (SSRSS_V2.cpp), SSRSStar_OptV2_Trsh (SSRSS_V3.cpp) functions. The output of these functions is a vector of size N, containing similarity of input query with all other nodes in the graph. Functions SSRSStar_OptV2 and SSRSStar_OptV2_Trsh, use an equivalency vector equiSet. Function FindEqueivalents computes the equiSet vector.

    /* ------ Finding Equivalents Nodes, IF USING SSRS* VERSION2 */
    vector<int> equiSet;
    FindEqueivalents(G, equiSet);

    /* ------ Single_Source Similarity computation ---------*/
    int query = 1;// query id
    vector<double> qSimVec(N, 0.2);//single source similarity vector
    SSRSStar_OptV1_Rec( G, indegVec, kmax, beta, lambda, query, qSimVec, count, theta, true);
    SSRSStar_OptV2_Rec( G, indegVec, kmax, beta, lambda, query, qSimVec, count, theta, true, equiSet);
    SSRSStar_OptV2_Trsh(G, indegVec, kmax, beta, lambda, query, qSimVec, count, theta, delta,  
                        true, equiSet, pruned_count);
