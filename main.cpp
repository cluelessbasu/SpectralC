#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "Triadic.h"
#include "Digraph.h"
#include "Graph.h"
#include "Digraph.h"
#include "GraphIO.h"
#include <set>
#include <list>

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
 std::list<std::set<VertexIdx>> output;
 std::vector<int> removals;
 

  if (loadGraph(argv[1], g, 1, IOFormat::escape))
    exit(1);

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
  std::vector<std::vector<VertexIdx>> adj(number, std::vector<VertexIdx>);
  std::vector<std::vector<VertexIdx>> hashmap(number, std::vector<VertexIdx>);

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

  while (counts){
    removals.clear();
    inf = cleaner(cg, degree, eps, adj);
    for (i=0; i<adj.size(); ++i){
        
        for (j=0; j<adj[i].size(); ++j){
            sources.push_back(i);
            dests.push_back(adj[i][j]);
            g_temp.srcs = &sources[0];
            g_temp.dsts = &dests[0];
        }

        if (adj[i].size()==0){
            adj.erase(adj.begin()+i);
            removals.push_back(i);
        }
    }

    int k = 0;

    for (i=0; i<removals.size(); ++i){
        if (!removals[k]>i)
            k++;
        else if (removals[k]!=i)
            hashmap[i-k].push_back(i);
    }

    g_temp.nVertices = adj.size();
    g_temp.nEdges = &sources[0].size();

    CGraph cg = makeCSR(g_temp);
    cg.sortById();
    CGraph cg_relabel = cg.renameByDegreeOrder();
    cg_relabel.sortById();
    CDAG dag = degreeOrdered(&cg_relabel);
    (dag.outlist).sortById();
    (dag.inlist).sortById();
    out = extractor(cg_relabel, degree, eps, inf, counts, adj);
    
    out_f.clear();

    for (i=0; i<out.size(); ++i ){
        out_f.insert(hashmap[i].back());
    }

    output.push_back(out_f);
  }


  FILE* f = fopen("out.txt","w");
  if (!f)
  {
      printf("could not write to output to out.txt\n");
      return 0;
  }
  fprintf(f,"%lld\n",cg_relabel.nVertices);
  fprintf(f,"%lld\n",cg_relabel.nEdges);
  for(int i = 0; i < 4; i++)
  {
      fprintf(f,"%f\n",nonInd[i]);
  }

  fclose(f);
}
