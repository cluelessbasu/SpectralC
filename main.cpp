#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "GraphIO.h"
#include "Graph.h"
#include "Triadic.h"
#include "Digraph.h"
#include "Digraph.h"
#include <set>
#include <list>
#include <cstdlib>
#include <cstdio>
//#include "FourVertex.h"
//#include "Conversion.h"
//#include "GetAllCounts.h"

using namespace Escape;
/*
void addEdge(vector<int> adj[], int u, int v)
{
    adj[u].push_back(v);
    adj[v].push_back(u);
}

// A utility function to delete an edge in an
// undirected graph.
void delEdge(vector<int> adj[], int u, int v)
{
    // Traversing through the first vector list
    // and removing the second element from it
    for (int i = 0; i < adj[u].size(); i++) {
        if (adj[u][i] == v) {
            adj[u].erase(adj[u].begin() + i);
            break;
        }  
    }
 
    // Traversing through the second vector list
    // and removing the first element from it
    for (int i = 0; i < adj[v].size(); i++) {
        if (adj[v][i] == u) {
            adj[v].erase(adj[v].begin() + i);
            break;
        }
    }
}
*/
int main(int argc, char *argv[])
{Graph g, g_temp;
 std::vector <VertexIdx> sources, dests;
 std::vector <VertexIdx> pVert;
 TriangleInfo2 inf;
 float eps = 0.1;
 VertexIdx i, j, k, l, v, vert, number, counts, *src, *dst;
 std::set<VertexIdx> out, out_f;
 std::list<std::set<VertexIdx> > output;
 std::vector<int> removals;
 long long int counter = 0;
 
 printf("Compiling started");

  if (loadGraph(argv[1], g, 1, IOFormat::escape)){
    printf("blah");
    exit(1);
  }
    

  printf("Loaded graph\n");
  CGraph cg = makeCSR(g);
  cg.sortById();
  printf("Converted to CSR\n");

  printf("Relabeling graph\n");
  CGraph cg_relabel = cg.renameByDegreeOrder();
  cg_relabel.sortById();
  printf("Creating DAG\n");
  CDAG dag = degreeOrdered(&cg_relabel);

  (dag.outlist).sortById();
  (dag.inlist).sortById();
  
  cg = cg_relabel.copy();

  number = cg.nVertices;
  std::vector<std::vector <VertexIdx> > adj(number);
  std::vector<std::vector <VertexIdx> > hashmap(number);

  double nonInd[4];
  src = g.srcs;
  dst = g.dsts;

  int degree[cg.nVertices];

  for (i = 0; i<sizeof(src)/sizeof(src[0]); ++i){
    j = src[i];
    k = dst[i];

    adj[j].push_back(k);
    adj[k].push_back(j);
  }
    
  for (i = 0; i<adj.size(); ++i){
    degree[i] = adj[i].size();
  }
  counts = cg.nVertices;
  FILE* f = fopen("out.txt","w");
  if (!f)
  {
      printf("could not write to output to out.txt\n");
      return 0;
  }
  while (counts){
    removals.clear();
    inf = cleaner(&cg, degree, eps, adj);
    for (i=0; i<adj.size(); ++i){
        
        for (j=0; j<adj[i].size(); ++j){
            sources.push_back(i);
            dests.push_back(adj[i][j]);
        }

        if (adj[i].size()==0){
            adj.erase(adj.begin()+i);
            removals.push_back(i);
        }
    }

    g_temp.srcs = &sources[0];
    g_temp.dsts = &dests[0];

    int k = 0;

    for (i=0; i<removals.size(); ++i){
        if (!removals[k]>i)
            k++;
        else if (removals[k]!=i)
            hashmap[i-k].push_back(i);
    }

    g_temp.nVertices = adj.size();
    g_temp.nEdges = sources.size();

    CGraph cg = makeCSR(g_temp);
    cg.sortById();
    CGraph cg_relabel = cg.renameByDegreeOrder();
    cg_relabel.sortById();
    CDAG dag = degreeOrdered(&cg_relabel);
    (dag.outlist).sortById();
    (dag.inlist).sortById();
    out = extractor(&cg_relabel, degree, eps, inf, counts, adj);
    
    out_f.clear();

    for (i=0; i<out.size(); ++i ){
        out_f.insert(hashmap[i].back());
        fprintf(f, "%lld,", hashmap[i].back());
    }
    counter++;
    output.push_back(out_f);
    fprintf(f, "End of cluster %lld,", counter);
  }



  fclose(f);
}
