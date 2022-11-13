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
#include <time.h>
#include <chrono>
#include <fstream>
 

//#include "FourVertex.h"
//#include "Conversion.h"
//#include "GetAllCounts.h"



using namespace Escape;

int main()
{Graph g, g_temp;

 auto t1 = time(NULL);
 std::vector <VertexIdx> sources, dests;
 std::vector <VertexIdx> pVert;
 TriangleInfo2 inf;
 float eps = 0.1;
 VertexIdx i, j, k, l, v, vert, number, counts, mapped;
 std::set<VertexIdx> out, out_f;
 //std::list<std::set<VertexIdx> > output;
 std::vector<int> removals, rem_temp;
 long int counter = 0, pos = 0;
 VertexIdx inter = 0, current = 0;
 VertexIdx avg_s=0, min_s=0, max_s=0;
 std::vector <VertexIdx> index;
 
 printf("Compiling started\n");
 loadGraph("Data/cit-Patents.edges", g, 1, IOFormat::escape);
  

  printf("Loaded graph\n");
  
  number = g.nVertices;

  printf("\n The number of vertices is %ld\n", number);

  std::vector<std::vector <VertexIdx> > adj(number);

  //printf("\nTesttttt");

  std::vector<std::vector <VertexIdx> > hashmap(number);
  
  printf("The length of the adjacency list is %ld", adj.size());
  fflush(stdin);
  double nonInd[4];

  VertexIdx* dst  = (VertexIdx*)malloc(sizeof(VertexIdx)*g.nEdges);
  VertexIdx* src  = (VertexIdx*)malloc(sizeof(VertexIdx)*g.nEdges);

  std::memcpy(src, g.srcs, sizeof(VertexIdx)*g.nEdges);
  std::memcpy(dst, g.dsts, sizeof(VertexIdx)*g.nEdges);

  
  int zeros = 0;
  
  //printf("\nShould repeat %ld times", sizeof(src)/sizeof(src[0]));
  int tester = 0;
  for (i = 0; i<g.nEdges; ++i){
    j = g.srcs[i];
    k = g.dsts[i];
    //sources.push_back(j);
    //dests.push_back(k);

    if (j>k)
        {
            adj[j].push_back(k);
            adj[k].push_back(j);}
    //printf("%ld, %ld\n", j,k);
  }

 delGraph(g);
 int* degree =  new int[number];
 for (i = 0; i<adj.size(); ++i){
    degree[i] = adj[i].size();
    if (degree[i]==0){
        //printf("%d",i);
        zeros++;
    }
  }
  index = sort_indexes (degree, adj.size());

 int edges = g.nEdges;
 std::cout << "\nMax is " << *std::max_element(degree, degree+g.nVertices)<< std::endl;
 std::cout << "Min is " << *std::min_element(degree, degree+g.nVertices)<< std::endl;


 printf("%d many zeros of %ld vertices\n",zeros, g.nVertices);

 counts =g.nVertices;
 std::ofstream files ("out.txt");
  
 std::vector <VertexIdx> flag (counts);
 std::vector <VertexIdx> reverse (counts);
 for (i=0; i<flag.size(); i++){
    flag[i] = 1;
    reverse[i] = i;
    hashmap[i].push_back(i);
 }

 //printf("\n i is %ld",i);

 int repeat = 0;
 printf("\n %ld vertices left, adj size: %ld ", counts, adj.size());
 auto t2 = time(NULL);
  printf("\n Clean begin at %s\n", ctime(&t2) );
  while (counts){
    repeat = 0;
   
    //printf("\n %ld vertices left", counts);
    while (repeat<6){
        
        repeat++;

        int remove_e = 0;
        rem_temp.clear();
        
        for (i=0; i<adj.size(); ++i){
        
                for (j=0; j<adj[i].size(); ++j){
                    sources.push_back(i);
                    //sources.push_back(adj[i][j]);
                    dests.push_back(adj[i][j]);
                    //dests.push_back(i);
                
            }
            /*      
            if (repeat > 1 && counter)
            {
                if (adj[i].size()==0){
                    //printf("\n i is %ld",i);
                    removals.push_back(hashmap[i].back());
                    rem_temp.push_back(i);
                    pos = hashmap[i].back();
                    flag[pos] = 0;
                    remove_e+= adj[i].size();
                    //printf ("\n Vetrex %ld removed", i);
            }
            }*/
        }
        if (sources.size()==0)
            break;
        clean(sources, dests, degree, eps, adj);
        
        //printf("\n%ld vertices  removed\n", rem_temp.size());
        
        /*
        for (i =0; i<rem_temp.size(); ++i){
            adj.erase(adj.begin()+rem_temp[i]-i);
            //printf("%ld ",removals[i]);
        }
        counts = adj.size();*/

        VertexIdx enums = 0;
        std::sort(removals.begin(), removals.end());
        std::sort(rem_temp.begin(), rem_temp.end());

        for (VertexIdx i =0; i<adj.size(); i++){
            enums += adj[i].size();
        }
        int k = 0;

        if (removals.size()>0){
            for (i=0; i<number; ++i){
                inter = removals[k];
                mapped = reverse[inter];
                current = reverse[i];

                if (mapped==i){

                    k++;
                    
                    }
                else {
                    if(mapped<i){
                        hashmap[i-k].pop_back();
                        hashmap[i-k].push_back(i);
                        reverse[i] = i-k;
                    }
                }
            }
        }
        
        
        

        //cg = makeCSR(g_temp);
        if (counter){
            repeat = 6;
            
        }

        sources.clear();
        dests.clear();
        rem_temp.clear();

        }
        t2 = time(NULL);
        printf ("\n Graph is cleaned, new graph has %ld vertices at %s\n", adj.size(), ctime(&t2));
        for (i=0; i<adj.size(); ++i){
        
                for (j=0; j<adj[i].size(); ++j){
                    sources.push_back(i);
                    //sources.push_back(adj[i][j]);
                    dests.push_back(adj[i][j]);
                    //dests.push_back(i);
                
            }}
            
/*
            if (counter || repeat > 1)
            {
                if (adj[i].size()==0){
                    //printf("\n i is %ld",i);
                    removals.push_back(hashmap[i].back());
                    rem_temp.push_back(i);
                    pos = hashmap[i].back();
                    flag[pos] = 0;
                    //printf ("\n Vetrex %ld removed", i);
            }
            }}
        
        for (i =0; i<rem_temp.size(); ++i){
            adj.erase(adj.begin()+rem_temp[i]-i);
            //printf("%ld ",removals[i]);
        }*/
        
        


        out = extract(adj, reverse,  flag, index);
        if (!counter)
            min_s = out.size();
        //printf ("\nWe got here: output size is %ld", out.size());
        if (out.size() && out.size()>max_s)
            max_s = out.size();
        else if (out.size() && out.size()<min_s)
            min_s = out.size();
        avg_s += out.size();
        out_f.clear();

        for (const VertexIdx &v :out){
            //out_f.insert(i);
            files<<v<<',';
            //printf("\n %ld,", hashmap[i].back());
            flag[v] = 0;
        }
        counter++;
        //counts -= out.size();
        //output.push_back(out_f);
        files<<"\nEnd of cluster "<< counter<<"\n";
        printf("\nEnd of cluster %ld with %ld vertices", counter, out.size()+1);
        //delete &g_temp, &cg, &cg_relabel;

        if (out.size()<2){
            VertexIdx edges = 0;
            for (i =0; i<adj.size(); i++){
                edges+= adj[i].size();
            }
            if (edges)
                break;
        }

        
    }

  printf("\n Min: %ld, Max %ld, Avg %ld, %ld clusters", min_s, max_s, avg_s/counter, counter);
  t2 = time(NULL);
  printf("\n%ld vertices  removed at %s\n", number-adj.size(),ctime(&t2));

  files.close();
  return 0;
}
