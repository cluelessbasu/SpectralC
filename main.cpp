//#include <algorithm>
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
#include <iterator>
#include <cstdint>
#include <cstring>
#include <numeric>
 

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

std::vector <VertexIdx> sort_indexes(const int* deg, long int size) {

  // initialize original index locations
  std::vector <VertexIdx> v(std::vector<VertexIdx>(deg, deg + size));
  std::vector<VertexIdx> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  
  

  return idx;
}



//int main(int argc, char *argv[])
int main()
{Graph g, g_temp;
 std::vector <VertexIdx> sources, dests;
 std::vector <VertexIdx> pVert;
 TriangleInfo2 inf;
 float eps = 0.1;
 VertexIdx i, j, k, l, v, vert, number, counts;
 std::set<VertexIdx> out, out_f;
 std::list<std::set<VertexIdx> > output;
 std::vector<int> removals;
 long int counter = 0, pos = 0;
 VertexIdx inter = 0;
 
 printf("Compiling started\n");
/*
  if (loadGraph(argv[1], g, 1, IOFormat::escape)){
    printf("blah");
    exit(1);
  }
*/
  loadGraph("ca-AstroPh.edges", g, 1, IOFormat::escape);
  

  printf("Loaded graph\n");
  CGraph cg = makeCSR(g);
  cg.sortById();
  printf("Converted to CSR\n");

  printf("Relabeling graph\n");
  CGraph cg_relabel; 
  cg_relabel  = cg.renameByDegreeOrder();
  //cg_relabel.sortById();
  //printf("Hehe\n");
  //printf("Creating DAG\n");
  //CDAG dag = degreeOrdered(&cg_relabel);

  //(dag.outlist).sortById();
  //(dag.inlist).sortById();
  
  cg = cg_relabel.copy();

  number = cg.nVertices;

  printf("\n The number of vertices is %ld\n", number);

  std::vector<std::vector <VertexIdx> > adj(number);

  //printf("\nTesttttt");

  std::vector<std::vector <VertexIdx> > hashmap(number);
  
  printf("The length of the adjacency list is %ld", adj.size());

  double nonInd[4];

  VertexIdx* dst  = (VertexIdx*)malloc(sizeof(VertexIdx)*g.nEdges);
  VertexIdx* src  = (VertexIdx*)malloc(sizeof(VertexIdx)*g.nEdges);

  std::memcpy(src, g.srcs, sizeof(VertexIdx)*g.nEdges);
  std::memcpy(dst, g.dsts, sizeof(VertexIdx)*g.nEdges);

  int degree[cg.nVertices];
  int zeros = 0;
  
  //printf("\nShould repeat %ld times", sizeof(src)/sizeof(src[0]));
  int tester = 0;
  for (i = 0; i<cg.nEdges; ++i){
    j = g.srcs[i];
    k = g.dsts[i];

    if (j>k)
        {
            adj[j].push_back(k);
            adj[k].push_back(j);}
    //printf("%ld, %ld\n", j,k);
  }
  printf("\n Number of edges is %ld \n", cg.nEdges);
  printf("\n Here's an example:\n");
  for (i = 0; i<adj[324].size(); ++i){
    printf ("%ld, ", adj[324][i]);
  }
  printf("\n");
  for (i = 0; i<adj[462].size(); ++i){
    printf ("%ld, ", adj[462][i]);
  }
  printf("\n");

  printf("\nTheir degrees are: %ld, %ld \n", adj[324].size(), adj[462].size());
  printf("\n Tester is %d\n",tester);

  for (i = 0; i<adj.size(); ++i){
    degree[i] = adj[i].size();
    if (degree[i]==0){
        //printf("%d",i);
        zeros++;
    }
  }

 std::cout << "Max is " << *std::max_element(degree, degree+cg.nVertices)<< std::endl;
 std::cout << "Min is " << *std::min_element(degree, degree+cg.nVertices)<< std::endl;


 std::vector <VertexIdx> indices = sort_indexes (degree, adj.size());

 printf("%d many zeros of %ld vertices\n",zeros, cg.nVertices);

  counts = cg.nVertices;
  FILE* f = fopen("out.txt","w");
  if (!f)
  {
      printf("could not write to output to out.txt\n");
      return 0;
  }

  
 std::vector <VertexIdx> flag (counts);
 std::vector <VertexIdx> reverse (counts);
 for (i=0; i<flag.size(); i++){
    flag[i] = 1;
    reverse[i] = i;
    hashmap[i].push_back(i);
 }
 
 //loop begins
  int repeat = 0;
  while (counts){
    repeat = 0;
   
    printf("\n %ld vertices left", counts);
    while (repeat<5){
        
        removals.clear();
        
        if (!counter)
            
            inf = cleaner(&cg, degree, eps, adj);
        else repeat = 5;

        repeat++;
        
        for (i=0; i<adj.size(); ++i){
        
            for (j=0; j<adj[i].size(); ++j){
                if (i<j){
                    sources.push_back(i);
                    dests.push_back(adj[i][j]);
                }
            }

            if (adj[i].size()==0){
                removals.push_back(i);
                pos = hashmap[i].back();
                flag[pos] = 0;
            }
        }

        printf("\n%ld vertices removed by cleaning\n", removals.size());
    
        for (i =0; i<removals.size(); ++i){
            adj.erase(adj.begin()+removals[i]-i);
            printf("%ld ",removals[i]);
        }
        counts = adj.size();

        //delete src, dst;

        VertexIdx* src = &sources[0]; 
        VertexIdx* dst = &dests[0];


        g_temp = newGraph(adj.size(), sources.size());

        //g_temp.nVertices = adj.size();
        //g_temp.nEdges = sources.size();

        //VertexIdx* dst  = (VertexIdx*)malloc(sizeof(VertexIdx)*g_temp.nEdges);
        //VertexIdx* src  = (VertexIdx*)malloc(sizeof(VertexIdx)*g_temp.nEdges);

        std::memcpy(g_temp.srcs, src, sizeof(VertexIdx)*g_temp.nEdges);
        std::memcpy(g_temp.dsts, dst, sizeof(VertexIdx)*g_temp.nEdges);

        int k = 0;

        for (i=0; i<number; ++i){
            inter = removals[k];
        

            if (hashmap[inter].back()<i){
                hashmap[i-k].pop_back();
                hashmap[i-k].push_back(i);
                reverse[i] = i-k;
                }
            else {
                if(hashmap[inter].back()==i){
                    k++;
                    //printf("\n %ld %ld %d %ld", inter, i, k, removals.size());
                }
            }
            }

        

        cg = makeCSR(g_temp);
        if (counter){
            inf = cleaner(&cg, degree, eps, adj);
        }
        }
        
        printf ("\n Graph is cleaned, new graph has %ld vertices\n", cg.nVertices);
        //cg.sortById();
        //cg_relabel = cg.renameByDegreeOrder();
        //cg_relabel.sortById();
        //dag = degreeOrdered(&cg_relabel);
        //(dag.outlist).sortById();
        //(dag.inlist).sortById();
        cg = cg_relabel.copy();
        out = extractor(&cg, degree, eps, inf, counts, adj, indices, flag, reverse, hashmap);
        printf ("\nWe got here: output size is %ld", out.size());
        out_f.clear();

        for (i=0; i<out.size(); ++i ){
            out_f.insert(hashmap[i].back());
            fprintf(f, "%ld,", hashmap[i].back());
            printf("\n %ld,", hashmap[i].back());
            flag[hashmap[i].back()] = 0;
        }
        counter++;
        //counts -= out.size();
        output.push_back(out_f);
        fprintf(f, "\nEnd of cluster %ld,", counter);
        printf("\nEnd of cluster %ld,", counter);
        //delete &g_temp, &cg, &cg_relabel;
    }



  fclose(f);
}
