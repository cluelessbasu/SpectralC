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
 VertexIdx i, j, k, l, v, vert, number, counts, mapped;
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
  //CGraph cg = makeCSR(g);
  //cg.sortById();
  //printf("Converted to CSR\n");

  //printf("Relabeling graph\n");
  //CGraph cg_relabel; 
  //cg_relabel  = cg.renameByDegreeOrder();
  //cg_relabel.sortById();
  //printf("Hehe\n");
  //printf("Creating DAG\n");
  //CDAG dag = degreeOrdered(&cg_relabel);

  //(dag.outlist).sortById();
  //(dag.inlist).sortById();
  
  //cg = cg_relabel.copy();

  number = g.nVertices;

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

  int degree[g.nVertices];
  int zeros = 0;
  
  //printf("\nShould repeat %ld times", sizeof(src)/sizeof(src[0]));
  int tester = 0;
  for (i = 0; i<g.nEdges; ++i){
    j = g.srcs[i];
    k = g.dsts[i];

    if (j>k)
        {
            adj[j].push_back(k);
            adj[k].push_back(j);}
    //printf("%ld, %ld\n", j,k);
  }

  int s1=0, s2=0;
        for (i =0; i<g.nEdges; i++){
            if (g.srcs[i]==0 && g.dsts[i]==0)
                s1++;
        }
  printf("\nHere are s1 and s2: %d, %d", s1, s2);
  printf("\n Number of edges is %ld \n", g.nEdges);
  printf("\n Here's an example:\n");
  
  for (i = 0; i<adj[324].size(); ++i){
    printf ("%ld, ", adj[324][i]);
  }
  printf("\n");
  printf("\n");
  for (i = 0; i<adj[462].size(); ++i){
    printf ("%ld, ", adj[462][i]);
  }
  printf("\n");

  printf("\nTheir degrees are: %ld, %ld \n", adj[324].size(), adj[462].size());
  
  for (i = 0; i<adj.size(); ++i){
    degree[i] = adj[i].size();
    if (degree[i]==0){
        //printf("%d",i);
        zeros++;
    }
  }

 printf("\nTheir degrees are: %ld, %ld \n", degree[324], degree[462]);
 int edges = g.nEdges;
 std::cout << "Max is " << *std::max_element(degree, degree+g.nVertices)<< std::endl;
 std::cout << "Min is " << *std::min_element(degree, degree+g.nVertices)<< std::endl;


 std::vector <VertexIdx> indices = sort_indexes (degree, adj.size());

 printf("%d many zeros of %ld vertices\n",zeros, g.nVertices);

  counts =g.nVertices;
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

 s1=0; 
 s2=0;
for (i =0; i<g.nEdges; i++){
    if (src[i]==0 && dst[i]==0)
        s1++;
}

printf("\nHere are s1 and s2: %d, %d", s1, s2);
 
 //loop begins

 //change to adj functions
  int repeat = 0;
  printf("\n Clean begin\n");
  while (counts){
    repeat = 0;
   
    printf("\n %ld vertices left", counts);
    while (repeat<5){
        
        removals.clear();
        
        if (!counter)
            
            clean(g, src, dst, degree, eps, adj);

        else repeat = 5;

        repeat++;

        int remove_e = 0;
        
        for (i=0; i<adj.size(); ++i){
        
            for (j=0; j<adj[i].size(); ++j){
                    sources.push_back(i);
                    sources.push_back(adj[i][j]);
                    dests.push_back(adj[i][j]);
                    dests.push_back(i);
                
            }

            if (adj[i].size()==0){
                removals.push_back(i);
                pos = hashmap[i].back();
                flag[pos] = 0;
                
            }
            remove_e+= adj[i].size();
        }

        printf("\n%ld vertices and %d edges removed by cleaning\n", removals.size(), edges - remove_e);

        edges-=(edges-remove_e);
    
        for (i =0; i<removals.size(); ++i){
            adj.erase(adj.begin()+removals[i]-i);
            //printf("%ld ",removals[i]);
        }
        counts = adj.size();

        //delete src, dst;

        VertexIdx* src = &sources[0]; 
        VertexIdx* dst = &dests[0];


        g = newGraph(adj.size(), sources.size());

        //g_temp.nVertices = adj.size();
        //g_temp.nEdges = sources.size();

        //VertexIdx* dst  = (VertexIdx*)malloc(sizeof(VertexIdx)*g_temp.nEdges);
        //VertexIdx* src  = (VertexIdx*)malloc(sizeof(VertexIdx)*g_temp.nEdges);

        std::memcpy(g.srcs, src, sizeof(VertexIdx)*g.nEdges);
        std::memcpy(g.dsts, dst, sizeof(VertexIdx)*g.nEdges);

        int k = 0;

        for (i=0; i<number; ++i){
            inter = removals[k];
            mapped = reverse[inter];

            if (mapped<i){

                hashmap[mapped-k].pop_back();
                hashmap[mapped-k].push_back(i);
                reverse[i] = mapped-k;
                }
            else {
                if(hashmap[mapped].back()==i){
                    k++;
                    //printf("\n %ld %ld %d %ld", inter, i, k, removals.size());
                }
            }
            }
        
        
        

        //cg = makeCSR(g_temp);
        if (counter){
            clean(g, src, dst, degree, eps, adj);
        }
        }
        
        printf ("\n Graph is cleaned, new graph has %ld vertices and %ld edges\n", g.nVertices, g.nEdges);
        //cg.sortById();
        //cg_relabel = cg.renameByDegreeOrder();
        //cg_relabel.sortById();
        //dag = degreeOrdered(&cg_relabel);
        //(dag.outlist).sortById();
        //(dag.inlist).sortById();
        //cg = cg_relabel.copy();
        out = extract(g.nVertices, adj, reverse, indices, flag);
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
