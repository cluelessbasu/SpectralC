#ifndef ESCAPE_TRIADIC_H_
#define ESCAPE_TRIADIC_H_

#include "ErrorCode.h"
#include "Graph.h"
#include "Digraph.h"
#include <list>
#include <set>
#include <vector>
#include <cstring>
#include <stdio.h>
using namespace Escape;

// Structure that stores triangle information of graph

struct TriangleInfo
{
    long int total; // total number of triangles
    long int *perVertex;  // array storing number of triangles for each vertex
    long int *perEdge;   // arry storing number of triangles for each edge

};

struct TriangleInfo2
{
    Count2 total; // total number of triangles
    Count2 *perVertex;  // array storing number of triangles for each vertex
    Count2 *perEdge;   // arry storing number of triangles for each edge

};


//use with a Graph or CGraph as argument
template <class T>
TriangleInfo newTriangleInfo(const T* graph)
{
  TriangleInfo ret;
  ret.total = 0;
  ret.perVertex = new Count[graph->nVertices];
  ret.perEdge = new Count[graph->nEdges];
  return ret;
}


void delTriangleInfo(TriangleInfo& info)
{
  delete[] info.perVertex;
  delete[] info.perEdge;
}


// Structure for storing all the triangles of graph
// Interpretation of the structure requires knowledge of the exact CGraph cg used to construct it
// Every index pos in the nbors list of cg corresponds to some edge (i,j). Note that each undirected edge appears as (i,j) and (j,i)
// Then trioffsets[pos] contains the starting index (in *triangles) of all the triangles that (i,j) is incident to.
// Basically, the portion of the array trioffsets[pos] to trioffsets[pos+1] in triangles is a list of vertices v1, v2,... such 
//that each vi forms a triangles with (i,j)
//
// We only store triangles for one copy of edge (i,j), where i < j in the degree ordering.


struct TriangleList
{
    VertexIdx total;
    VertexIdx *triangles;
    EdgeIdx *trioffsets;
};

std::set<long int> Union(const std::set<long int> & s0, const std::set<long int> & s1) {
    std::set<long int> s; 
    std::set_union(s0.begin(), s0.end(),
        s1.begin(), s1.end(),
        std::inserter(s,s.begin()));

    return s;
}



// The wedge enumeration algorithm that produces all triangles
// Input: a pointer g to a CGraph, that CAN be a DAG. Indeed, we will call wedgeEnumerator on DAGs.
// Output: a TriangleInfo for g. The ordering of edges in perEdge (of TriangleInfo) is that same as g.

TriangleInfo wedgeEnumerator(CGraph *g)
{
   TriangleInfo ret;   // output 
   ret.total = 0;      // initialize outout
   ret.perVertex = new Count[g->nVertices+1];
   ret.perEdge = new Count[g->nEdges+1]; 

   for (VertexIdx i=0; i < g->nVertices; ++i)
       ret.perVertex[i] = 0;

   for (EdgeIdx j=0; j < g->nEdges; ++j)
       ret.perEdge[j] = 0; 

   VertexIdx end1, end2;
   EdgeIdx loc1, loc2;

   for (VertexIdx i=0; i < g->nVertices; ++i) // loop over vertices
       for (EdgeIdx j = g->offsets[i]; j < g->offsets[i+1]; ++j)   // loop over neighbor of i
           for (EdgeIdx k = j+1; k < g->offsets[i+1]; ++k)         // loop over another neighbor of i
           {
               end1 = g->nbors[j];     // we are now looking at wedge (i, end1, end2), centered at i
               end2 = g->nbors[k];

               loc1 = g->isEdge(end1,end2);  // check if edge (end1, end2) is present
               loc2 = g->isEdge(end2,end1);  // check if (end2, end1) is present. Could be either, since g is not necessarily undirected
               if (loc1 != -1)        // (end1, end2) is present
               {
                   ret.total++;       // found a triangle! So update total.

                   ret.perVertex[i]++; // update all per vertex counts
                   ret.perVertex[end1]++;
                   ret.perVertex[end2]++;

                   ret.perEdge[j]++;  // update all per edge counts. Note that location used is same as position in g->nbors
                   ret.perEdge[k]++;
                   ret.perEdge[loc1]++;
               }
               if (loc2 != -1)       // (end2, end1) is present
               {
                   ret.total++;      // found a triangle! Update total.

                   ret.perVertex[i]++; // update all per vertex counts
                   ret.perVertex[end1]++;
                   ret.perVertex[end2]++;

                   ret.perEdge[j]++;  // update all per edge counts
                   ret.perEdge[k]++;
                   ret.perEdge[loc2]++;
               }


           }

   return ret;
}


// This wedge enumeration algorithm produces all triangles, and is more
// efficient than the previous version. We use binary search to find edges, and also pass
// the original graph to cut down on queries.
// Input: a pointer gout to a CGraph labeled according to degree
// Output: a TriangleInfo for g. The ordering of edges in perEdge (of TriangleInfo) is that same as g.

void clean (std::vector <VertexIdx>& src, std::vector <VertexIdx>& dst, int *degree, float eps, std::vector<std::vector<VertexIdx> >& adj){
        //VertexIdx enums = G.nEdges;
        //VertexIdx vnum = G.nVertices;
        float twt = 0;
        VertexIdx u, v;

        VertexIdx vnum = adj.size();
        VertexIdx enums = src.size();

        //VertexIdx sources[enums], dests[enums];
        float weights[enums];

        for (VertexIdx i = 0; i<enums; i++){
            twt = 0;
            u = src[i];
            v = dst[i];
            //printf("\n Edges is %ld %ld \n", u,v);
            for (VertexIdx j =0; j<adj[u].size(); j++){
                for (VertexIdx k =0; k<adj[v].size(); k++){
                    if (adj[u][j] == adj[v][k]){
                        twt += 1.0/degree[adj[v][k]];
                        //printf("\nTriwt is %f, degree is %d ",twt, degree[adj[v][k]]);
                    }
                }
            }
            if (twt<eps){
                //printf("Flag\n");
                for (VertexIdx j =0; j<adj[u].size(); j++){
                    if (adj[u][j] == v)
                        adj[u].erase(adj[u].begin()+j);
                }
                for (VertexIdx j =0; j<adj[v].size(); j++){
                    if (adj[v][j] == u)
                        adj[v].erase(adj[v].begin()+j);
                        
                }
            }
        }
}

std::set <VertexIdx> extract (std::vector<std::vector<VertexIdx> >& adj, std::vector <VertexIdx>& rev, std::vector <VertexIdx>& indices, std::vector <VertexIdx>& flag){

        VertexIdx i, v1, v2, v, u; 
        VertexIdx vnum = adj.size();
        std::set <VertexIdx> set1;
        for (VertexIdx blah = 0; blah < vnum; blah++){
            while(!(flag[indices[blah]]) || adj[rev[indices[blah]]].size()<3)
                            {blah++;}
            VertexIdx i_n = indices[blah];

            i = rev[i_n];
            set1.insert(i);

            for (VertexIdx j= 0; j<adj[i].size(); j++){
                if (adj[j].size()<2*adj[i].size())
                    set1.insert(j);
            }

            for (EdgeIdx j = 0; j<set1.size(); ++j){
                v1 = *std::next(set1.begin(), j);
                for (EdgeIdx k = 0; k<adj[v1].size(); k++){
                    set1.insert(adj[v1][k]);
                }
            }

            for (EdgeIdx j = 0; j<set1.size(); ++j){
                int deg = 0;
                v1 = *std::next(set1.begin(), j);
                for (EdgeIdx k = 0; k<set1.size(); ++k){
                    v2 = *std::next(set1.begin(), k);
                    for (EdgeIdx l=0; l<adj[v1].size(); ++l){
                        if (adj[v1][l] == v2)
                            deg++;
                    }
                    if (!(deg>adj[v1].size()/2)){
                        set1.erase(v1);
                    }
                }

            }

        for (VertexIdx i = 0; i < adj.size(); i++){

               for (EdgeIdx j = 0; j<set1.size(); j++){

                    v = *std::next(set1.begin(), j);

                    for (VertexIdx k = 0; k<adj[i].size(); k++){
                        if (v == adj[i][k])
                            adj[i].erase(adj[i].begin()+k);
                    }

               }

        } 


        for (EdgeIdx j = 0; j<set1.size(); j++){

                v = *std::next(set1.begin(), j);
                adj.erase(adj.begin()+v);
                flag[rev[v]] = 0;

        }

        
    return set1;
    }    
}

TriangleInfo2 cleaner(CGraph *gout, int *degree, float eps, std::vector<std::vector<VertexIdx> >& adj) //fix output format
{  std::set <VertexIdx> set1;
   TriangleInfo2 ret;   // output 
   ret.total = 0;      // initialize outout
   ret.perVertex = new Count2[gout->nVertices+1];
   ret.perEdge = new Count2[gout->nEdges+1]; 
   float weight;

   for (VertexIdx i=0; i < gout->nVertices; ++i)
       ret.perVertex[i] = 0;

   for (EdgeIdx j=0; j < gout->nEdges; ++j)
       ret.perEdge[j] = 0; 

   VertexIdx end1, end2;
   EdgeIdx loc;

   for (VertexIdx i=0; i < gout->nVertices; ++i) // loop over vertices
       for (EdgeIdx j = gout->offsets[i]; j < gout->offsets[i+1]; ++j)   // loop over neighbor of i
           for (EdgeIdx k = j+1; k < gout->offsets[i+1]; ++k)         // loop over another neighbor of i
           {
               end1 = gout->nbors[j];     // we are now looking at wedge (i, end1, end2), centered at i
               end2 = gout->nbors[k];

               // note that end1 < end2 because of the labeled ordering

               loc = gout->getEdgeBinary(end1,end2);
               if(degree[i]*degree[end1]*degree[end2]!=0){
                    if (loc != -1)        // (end1, end2) is present
                    {   
                        ret.total++;       // found a triangle! So update total.
                        weight = 1/(degree[i]*degree[end1]*degree[end2]);

                        ret.perVertex[i] += weight; // update all per vertex counts
                        ret.perVertex[end1] += weight;
                        ret.perVertex[end2] += weight;

                        ret.perEdge[j] += weight;  // update all per edge counts. Note that location used is same as position in g->nbors
                        ret.perEdge[k] += weight;
                        ret.perEdge[loc] += weight;
                    }

                    if (ret.perEdge[k]< eps/(degree[end1]*degree[end2])){
                        //edgeFlag[k] = 0;
                                adj[end1].erase(std::remove(adj[end1].begin(), adj[end1].end(), end2), adj[end1].end());
                                adj[end2].erase(std::remove(adj[end2].begin(), adj[end2].end(), end1), adj[end2].end());
                                //if len(adj[end1]==0){
                                //    vtxFlag[end1] = 0;
                                //}
                                //if len(adj[end2]==0){
                                //     vtxFlag[end2] = 0;
                                //}
                        
               }
               }

               
           }
    
   

   return ret;
}

std::set <VertexIdx> extractor(CGraph *gout, int *degree, float eps, TriangleInfo2 info, VertexIdx &count, std::vector<std::vector<VertexIdx> >& adj, std::vector <VertexIdx>& indices, std::vector <VertexIdx>& flag, std::vector <VertexIdx>& rev, std::vector<std::vector<VertexIdx> >& hash){

        VertexIdx end1, end2, e1, a,b,c,d, i, j, l, v, v_r, k, triangle;
        EdgeIdx loc, loc1, loc2, loc3;

        std::set <VertexIdx> set1;
        std::set <VertexIdx> set2;
        std::set <VertexIdx> out;

        float incident = 0;
        float inside = 0;
        float weight = 0;
        float wt = 0;

        k = 0;

        TriangleInfo2 ret;  
        

        for (VertexIdx blah=0; blah < flag.size(); ++blah) // loop over vertices
            {           
                        while(!(flag[indices[blah]]) || adj[rev[indices[blah]]].size()<3)
                            {blah++;}
                        VertexIdx i_n = indices[blah];

                        i = rev[i_n];
                        printf("\n%ld, %ld", (flag[indices[blah]]) , adj[i].size());
                        printf ("\nVertex is %ld/%ld, degree is %d, blah is %ld\n ", i_n,i,degree[i_n], blah);

                        //VertexIdx next = 0;
                        //blah++;
                        //while(!(flag[indices[blah]]) || adj[indices[blah]].size()<3)
                        //    {blah++;}

                        //VertexIdx next_n = indices[blah];
                        //next = rev[next_n];

                        //printf("The next one is %ld with degree %d, the offsets are %ld and %ld\n ", next, degree[next], gout->offsets[i], gout->offsets[next]);

                        set1.insert(i);

                        printf ("Offsets are %ld, %ld",gout->offsets[i+1], gout->offsets[i] );
                
                        for (EdgeIdx j = gout->offsets[i]; j < gout->offsets[i+1]; ++j)   // loop over neighbor of i
                        {   
                            v = gout->nbors[j];
                            v_r = rev[v];
                            printf ("\nDegree of this nbr is %d\n", degree[v_r]);
                            printf ("\n Compare %ld and %ld", gout->offsets[i+1]-gout->offsets[i], adj[i].size());
                            if (degree[v_r]<degree[i_n]/(pow(eps,3))){
                                //k ++;
                                set1.insert(v);
                               
                            }
                            else {
                                printf ("\n Degree of %ld is %d\n", v_r, degree[v_r]);
                                printf ("\n Compare %ld and %ld", gout->offsets[i+1]-gout->offsets[i], adj[v_r].size());
                            }
                        }//printf("\nExtractor loop 1");
                        printf("\n Set 1 has %ld elements\n", set1.size());
                        for (EdgeIdx j = 0; j<set1.size(); ++j){
                            e1 = *std::next(set1.begin(), j);
                            //e1 = set1[j];
                            a = gout->nbors[e1];
                            set2.insert(a);

                        }//printf("\nExtractor loop 2");
                        ret.total = 0;      // initialize outout
                        ret.perVertex = new Count2[set2.size()];
                        ret.perEdge = new Count2[set2.size()];
                        VertexIdx end1, end2;
                        EdgeIdx loc;

                        for (VertexIdx ii=0; ii < set2.size(); ++ii) // loop over vertices
                            {
                                i = *std::next(set2.begin(), ii);
                                //i = set2[ii];
                                for (EdgeIdx j = gout->offsets[i]; j < gout->offsets[i+1]; ++j)   // loop over neighbor of i
                                    for (EdgeIdx k = j+1; k < gout->offsets[i+1]; ++k)         // loop over another neighbor of i
                                    {
                                        end1 = gout->nbors[j];     // we are now looking at wedge (i, end1, end2), centered at i
                                        end2 = gout->nbors[k];

                                        // note that end1 < end2 because of the labeled ordering

                                        loc = gout->getEdgeBinary(end1,end2);
                                        //triangle = gout->trioffsets[k];
                                        if (loc != -1 && set2.count(i) && set2.count(end1) && set2.count(end2) && loc)        // (end1, end2) is present
                                        {   
                                            wt = 1/(degree[rev[i]]*degree[rev[end1]]*degree[rev[end2]]);

                                            ret.perVertex[ii] += weight; // update all per vertex counts
                                            ret.perVertex[end1] += weight;
                                            ret.perVertex[end2] += weight;
                                        }

                                        }
                                    }//printf("\nExtractor loop 3");

                            long int values[set2.size()][3];



                            for (i = 0 ; i < set2.size(); ++i){
                                values[i][0] = *std::next(set2.begin(), i);
                                //values[i][0] = set2[i];
                                values[i][1] = i;
                                values[i][2] = ret.perVertex[i];
                                incident +=  info.perVertex[*std::next(set2.begin(), i)];
                            }//printf("\nExtractor loop 4");

                            for (i=0; i<set2.size(); ++i){
                                for ( VertexIdx j = i; j > 0; j--)
                                    if (values[j][2] > values[j-1][2]) {
                                        std::swap(values[j][0], values[j-1][0]);
                                        std::swap(values[j][1], values[j-1][1]);
                                        std::swap(values[j][2], values[j-1][2]);

                                    }
                            }//printf("\nExtractor loop 5");

                            for (i = 0; i <set2.size(); ++i){
                                out.insert(values[i][0]);
                                inside += ret.perVertex[values[i][1]];
                                if (inside >= (pow(eps,2))*incident){
                                    break;
                                }

                            }//printf("\nExtractor loop 6");

                            out = Union(set1, out);
                            printf("\n Right now out is %ld", out.size());

                            for (VertexIdx i=0; i < gout->nVertices; ++i) // loop over vertices
                                for (EdgeIdx j = gout->offsets[i]; j < gout->offsets[i+1]; ++j)   // loop over neighbor of i
                                    for (EdgeIdx k = j+1; k < gout->offsets[i+1]; ++k)         // loop over another neighbor of i
                                    {
                                        end1 = gout->nbors[j];     // we are now looking at wedge (i, end1, end2), centered at i
                                        end2 = gout->nbors[k];

                                        // note that end1 < end2 because of the labeled ordering

                                        a = std::min(i, end1);
                                        b = std::max(i, end1);
                                        c = std::min(i, end2);
                                        d = std::min(i, end2);

                                        loc1 = gout->getEdgeBinary(end1,end2);
                                        loc2 = gout->getEdgeBinary(a,b);
                                        loc3 = gout->getEdgeBinary(c,d);
                                        /*
                                        if (loc1 != -1 && (out.find(end1)! = out.end() || out.find(end2)! = out.end()) )        
                                            edgeFlag[loc1] = 0;
                                        if (loc2 != -1 && (out.find(a)! = out.end() || out.find(b)! = out.end()) )        
                                            edgeFlag[loc2] = 0;
                                        if (loc3 != -1 && (out.find(c)! = out.end() || out.find(d)! = out.end()) )        
                                            edgeFlag[loc3] = 0;
                                        */

                                        

                                        
                                    }//printf("\nExtractor loop 7");

                            int dummy =0;

                            for (i = 0; i<adj.size(); ++i){
                                
                                for (long int vert = 0; vert<out.size(); ++vert){
                                    
                                    v = *std::next(out.begin(), vert);
                                    //printf("\n To remove: %ld", v);
                                    if (std::find(adj[i].begin(), adj[i].end(), v) != adj[i].end()){
                                        //adj[i].erase(std::remove(adj[i].begin(), adj[i].end(), v), adj[i].end());
                                        std::vector <VertexIdx>::iterator dels; 
                                         dels = std::remove(adj[i].begin(),adj[i].end(),v);
                                         //printf("\nFOUND: %ld", v);
                                        
                                        dummy++;
                                    }
                                }
                            }//printf("\nExtractor loop 8");
                        break;
                        }//printf("\nExtractor int loop");
                        
                   //break;      
        
             //printf("\nExtractor out loop");
        printf("\nEnd of extractor\n");
        count = count - out.size();
        return out;  // possibly return the net incident triangle weight to reduce work later on: what about the edge list? 

}


// 
// ccPerDeg: This function computes the degree-wise clustering coefficients. 
// Input: 
//     g: A pointer to the original graph
//     gout: A pointer to a Degree-Ordered DAG, built from g, that is used for the actual triangle counting
//     ccdegarray: A pointer to array of floats that contains the output. This must be of length at least the number of vertices.
//      The entry ccperdegarray[deg] is the average clustering coefficient of degree deg vertices.
//      If there are no vertices of a given degree, it is vacuously set to zero.
// Output:
//     Count: the maximum degree (so one only needs to read that many entries of ccarray
// 
// The function works by first calling betterWedgeEnumerator on gout, to get the TriangleInfo structure
// with per vertex counts. Then, the function computes the clustering coefficient for each vertex.
// The degree of the vertex is obtained from g (which is why it is passed). The clustering coefficients
// are appropriately binned to get the final output.
//     

Count ccPerDeg(CGraph *g, CGraph *gout, float *ccdegarray)
{
    VertexIdx maxdeg = 0; //for maximum degree
    VertexIdx i; //indices for looping
    VertexIdx deg; //for storing degree
    VertexIdx *degdist; // array for storing degree distribution
    TriangleInfo info; //for storing the triangle info

    degdist = new VertexIdx[g->nVertices+1]; // initializing array for degree distribution
    info = wedgeEnumerator(gout); // get all the triangle info, from wedge enumeration on the DOG DAG gout
    
    for (i=0; i<g->nVertices; ++i) //loop over vertices
    {
        degdist[i] = 0; //initialize to zero
        ccdegarray[i] = 0.0; //initialize to zero
    }

    for (i=0; i<g->nVertices; ++i) //loop over vertices
    {
        deg = g->offsets[i+1] - g->offsets[i]; // compute degree of vertex i
        if (deg > maxdeg)
            maxdeg = deg; // update maximum degree
//         printf("%lld %lld %f\n",info.perVertex[i],deg,info.perVertex[i]*2/(float)(deg*(deg-1)));
        if (deg > 1) // only do for deg > 1, otherwise there is a divide by zero
            ccdegarray[deg] += (float) info.perVertex[i]*2/(float)(deg*(deg-1)); //adding clustering coefficient of vertex, at appropriate position in ccdegarray
        degdist[deg]++; //updating number of vertices of degree deg
    }

    // at this point, ccperdegarray has the *sum* of clustering coefficients, per degree. degdist has the degree counts,
    // so it only remains to divide to get the final ccperdegarray
    for (deg=1; deg<=maxdeg; ++deg) // loop over degrees up to maxdeg
    {
        if (degdist[deg] == 0) // if no vertices of this degree, continue
            continue;
        ccdegarray[deg] = ccdegarray[deg]/(float) degdist[deg]; // compute the average clustering coefficient per degree
    }

    delete[] degdist; // free memory used 
    delTriangleInfo(info); // free triangle info

    return maxdeg; // return maximum degree, as promised
}

// This assumes that gout and gin are *reverses* of each other
// Given the triangle info with respect to a DAG gout, this outputs the triangle info with respect to the reversed DAG gin.

TriangleInfo moveOutToIn(CGraph *gout, CGraph *gin, TriangleInfo *outinfo)
{
    TriangleInfo ret;
    VertexIdx i,j;
    EdgeIdx loc;
    ret.perVertex = new EdgeIdx[gout->nVertices+1];
    ret.perEdge = new EdgeIdx[gout->nEdges+1]; 
 
    ret.total = outinfo -> total; // total triangles the same 

    for (i=0; i < gout->nVertices; ++i)
    {
        ret.perVertex[i] = outinfo->perVertex[i]; // vertex triangle count the same
        for (EdgeIdx pos = gout->offsets[i]; pos < gout->offsets[i+1]; ++pos)
        {
            j = gout->nbors[pos]; // look at edge i->j
            loc = gin->getEdgeBinary(j,i);

            if (loc == -1) // Edge (j,i) not in gin, so arguments are not reverses of each other
            {
                printf("Error in moveOutToIn: gout and gin not reverses of each other\n");
                printf("i = %ld, j = %ld\n",i,j);
                printf("%ld: ",i);
                for (EdgeIdx posnew = gout->offsets[i]; posnew < gout->offsets[i+1]; posnew++)
                    printf("%ld ",gout->nbors[posnew]);
                printf("\n");
                printf("%ld: ",j);
                for (EdgeIdx posnew = gin->offsets[j]; posnew < gin->offsets[j+1]; posnew++)
                    printf("%ld ",gin->nbors[posnew]);
                printf("\n");
                exit(EXIT_FAILURE);
            }
            ret.perEdge[loc] = outinfo->perEdge[pos]; // loc is index in gin for edge (j,i), so we assign the value of triangle count for (i,j)
        }
    }

    return ret;
}


// This stores all triangles in a TriangleList structure, corresponding to CGraph g. 
// The number of triangles numtri is required for initialization

TriangleList storeAllTriangles(CGraph *g, EdgeIdx numtri)
{
    TriangleList ret;

    ret.total = 3*numtri;
    ret.triangles = new VertexIdx[3*numtri+1];
    ret.trioffsets = new EdgeIdx[g->nEdges+1];

    EdgeIdx posj = 0;

    EdgeIdx current = 0;
    for (VertexIdx i=0; i < g->nVertices; i++)
    {
        VertexIdx degi = g->offsets[i+1] - g->offsets[i];
        for (posj = g->offsets[i]; posj < g->offsets[i+1]; posj++)
        {
            VertexIdx j = g->nbors[posj];
//             if (i == g->nVertices-1)
//                 printf("j is %lld\n",j);
            VertexIdx degj = g->offsets[j+1] - g->offsets[j];
            ret.trioffsets[posj] = current;
            if (degj < degi || (degj == degi && j <= i))
                continue;
            
            for (EdgeIdx ptr = g->offsets[i]; ptr < g->offsets[i+1]; ptr++)
            {
                VertexIdx nbr = g->nbors[ptr];
                if (g->getEdgeBinary(nbr,j) != -1)
                {
//                     if (i == g->nVertices-1)
//                         printf("%lld at ptr %lld: nbr is %lld\n",current,ptr,nbr);
                    ret.triangles[current] = nbr;
                    current++;
                }
            }
        }
    }
    
    ret.trioffsets[posj] = current; // posj is the index after the last edge in g. we set this final offset to the current position

    return ret;
}

// c-triangle-pruning is achieved by removingevery edge that participates in less than c triangles
// Input: CGraph outDAG gout, a value c
// Output: CGraph representation of the c-triangle-pruned graph

CGraph cTriPrune(CGraph *gout, VertexIdx c)
{
    CGraph ret;    // final return value
    Graph retEdges;  // initially, we'll construct truss as list of edges

    retEdges.nVertices = gout->nVertices;
    retEdges.nEdges = 0;
    retEdges.srcs = new VertexIdx[2*gout->nEdges+1];
    retEdges.dsts = new VertexIdx[2*gout->nEdges+1];
    TriangleInfo info = wedgeEnumerator(gout);

    for (VertexIdx i=0; i < gout->nVertices; i++)
        for (EdgeIdx posj=gout->offsets[i]; posj < gout->offsets[i+1]; posj++) //looping over all edges
        {
            VertexIdx j = gout->nbors[posj]; // edge (i,j)
            if(info.perEdge[posj] >= c) // edge (i,j) has more than c triangles
            {
               retEdges.srcs[retEdges.nEdges] = i; // insert edge (i,j) in retEdges
               retEdges.dsts[retEdges.nEdges] = j;
               retEdges.nEdges++;

               retEdges.srcs[retEdges.nEdges] = j; // insert edge (j,i) in retEdges
               retEdges.dsts[retEdges.nEdges] = i;
               retEdges.nEdges++;
            }
        }

   ret = makeCSR(retEdges);
   return ret;
}


// A c-truss is the largest subgraph where every edge participates in at least c triangles
// Input: CGraph g, the associated TriangleList tlist, a value c
// Output: CGraph representation of the c-truss

CGraph cTruss(CGraph *g, TriangleList* tlist, VertexIdx c)
{
    CGraph ret;    // final return value
    Graph retEdges;  // initially, we'll construct truss as list of edges

    EdgeIdx *triCount = new EdgeIdx[g->nEdges+1];  // store triangle count of each edge
    EdgeInfo *toDelete = new EdgeInfo[g->nEdges+1]; // store list of edges to be deleted

    EdgeIdx ind_toDelete = 0; // largest index in toDelete

    bool *deleted = new bool[g->nEdges+1];  // store flags for deleted edges

    for (EdgeIdx i=0; i < g->nEdges+1; i++) // initialize all edges as not deleted
        deleted[i] = false;

    for (VertexIdx i=0; i < g->nVertices; i++)
        for (EdgeIdx posj=g->offsets[i]; posj < g->offsets[i+1]; posj++) //looping over all edges
        {
            VertexIdx j = g->nbors[posj]; // edge (i,j)
            if (i > j)   // only look at pairs ordered properly
                continue;
            triCount[posj] = tlist->trioffsets[posj+1] - tlist->trioffsets[posj];  // store the number of triangles that edge (i,j) participates in. Information is exactly stored in tlist
//             printf("%lld %lld %lld\n",i,j,triCount[posj]);
            if (triCount[posj] < c) // edge should be deleted
            {
                EdgeInfo current;
                current.src = i;
                current.dest = j;
                current.index = posj;
                toDelete[ind_toDelete] = current;   // storing current edge in toDelete
                ind_toDelete++;  // update index
            }
        }

    while(ind_toDelete >= 1) // while toDelete is non-empty
    {
        EdgeInfo toRemove = toDelete[ind_toDelete-1]; // get edge to remove
        ind_toDelete--; // update last index in toDelete
        VertexIdx i = toRemove.src;
        VertexIdx j = toRemove.dest;
        deleted[toRemove.index] = true;
       
//         printf("%lld %lld\n",i,j); 
        for (EdgeIdx indk = tlist->trioffsets[toRemove.index]; indk < tlist->trioffsets[toRemove.index+1]; indk++) // looping over triangles that edge participates in
        {
            VertexIdx k = tlist->triangles[indk]; // (i,j,k) forms triangle

            EdgeIdx locik, locjk;
            if (i < k) // get position of ordered edge (i,k)
                locik = g->getEdgeBinary(i,k);
            else
                locik = g->getEdgeBinary(k,i);

            if (locik == -1) // something is wrong. (i,k) has to be edge in g
            {
                printf("Error: edge i,k is not present in graph, but triangle (i,j,k) is stored in tlist\n");
                exit(1);
            }
            if (deleted[locik])  // if (i,k) has been already deleted, then this triangle has been deleted, so continue
                continue;

            if (j < k) // get position of ordered edge (j,k)
                locjk = g->getEdgeBinary(j,k);
            else
                locjk = g->getEdgeBinary(k,j);

            if (locjk == -1) // something is wrong. (j,k) must be edge in g
            {
                printf("Error: edge j,k not present in graph, but triangle (i,j,k) is stored in tlist\n");
                exit(1);
            }
            if (deleted[locjk]) // if (j,k) has been deleted, triangle (i,j,k) has been deleted, so continue
                continue;


            // we delete triangle (i,j,k), so decrement triangle counts appropriately
            triCount[locik]--; // decrement count for (i,k)
            if (triCount[locik] < c) // (i,k) now participates in less than c triangles
            {
                EdgeInfo next = {i,k,locik}; // we should delete (i,k)
                toDelete[ind_toDelete] = next;
                ind_toDelete++;
            }
            triCount[locjk]--; // decrement count for (j,k)
            if (triCount[locjk] < c) // (j,k) now participates in less than c triangles
            {
                EdgeInfo next = {j,k,locjk}; // we should delete (j,k)
                toDelete[ind_toDelete] = next;
                ind_toDelete++;
            }
        }
    }
   // now construct graph with non-deleted edges

   retEdges.nVertices = g->nVertices;
   retEdges.nEdges = 0;
   retEdges.srcs = new VertexIdx[g->nEdges+1];
   retEdges.dsts = new VertexIdx[g->nEdges+1];
   for (VertexIdx i = 0; i < g->nVertices; i++)
       for (EdgeIdx posj = g->offsets[i]; posj < g->offsets[i+1]; posj++) // loop through edges in g
       {
           VertexIdx j = g->nbors[posj];
           if (i < j && !deleted[posj]) // if (i,j) is ordered and not deleted
           {
               retEdges.srcs[retEdges.nEdges] = i; // insert edge (i,j) in retEdges
               retEdges.dsts[retEdges.nEdges] = j;
               retEdges.nEdges++;

               retEdges.srcs[retEdges.nEdges] = j; // insert edge (j,i) in retEdges
               retEdges.dsts[retEdges.nEdges] = i;
               retEdges.nEdges++;
           }
       }

   ret = makeCSR(retEdges);
   return ret;
}

// This function computes closure rate as a function of common neighbors. Consider all pairs (i,j) that
//have exactly c neighbors in common. The closure rate for c is the fraction of such pairs that are also edges.
// Input: CGraph g, array common, array closed
// Output: Final length of arrays common and closed. These arrays are populated with desired output. ith element of common is the number of pairs of vertices that have i neighbors in common.
// The ith element of closed is the number of such pairs that are also edges. 

Count cClosure(CGraph* g, Count* common, Count* closed)
{
    Count ret = 0; // this will eventually be the maximum number of common neighbors
    VertexIdx *wedge_count = new VertexIdx[g->nVertices+1];

    for (VertexIdx i=0; i < g->nVertices; i++) // initialize arrays to 0
    {
        common[i] = 0;
        closed[i] = 0;
        wedge_count[i] = 0;
    }

    for (VertexIdx i=0; i < g->nVertices; i++) // loop over vertices
    {
        for (EdgeIdx posj = g->offsets[i]; posj < g->offsets[i+1]; posj++) // loop over neighbors of i
        {
            VertexIdx j = g->nbors[posj];
            for (EdgeIdx posk = g->offsets[j]; posk < g->offsets[j+1]; posk++)
            {
                VertexIdx k = g->nbors[posk]; // i-j-k is wedge
                if (k <= i) //k is lower in order, so ignore to prevent double counting
                    continue;
                wedge_count[k]++; // update number of wedges ending at k
            }
        }

        for (EdgeIdx posj = g->offsets[i]; posj < g->offsets[i+1]; posj++) // loop over neighbors of i
        {
            VertexIdx j = g->nbors[posj];
            for (EdgeIdx posk = g->offsets[j]; posk < g->offsets[j+1]; posk++)
            {
                VertexIdx k = g->nbors[posk]; // i-j-k is wedge
                if (k <= i) //k is lower in order, so ignore to prevent double counting
                    continue;
                if (wedge_count[k] != 0) // we have not processed (i,k)
                {
                    common[wedge_count[k]]++; // there are exactly wedge_count[k] common neighbors between i and k
                    if (g->isEdgeBinary(i,k)) // (i,k) is edge, so this pair is closed
                        closed[wedge_count[k]]++; // update closed count
                    if (wedge_count[k] > ret)
                        ret = wedge_count[k]; // update the maximum number of common neighbors
                    wedge_count[k] = 0; // reset
                }
            }
        }
    }

    return ret;
}

// Debug function. Prints out data in a TriangleInfo structure. Uses CGraph g to enumerate properly.

void printTri(FILE* f, TriangleInfo *info, CGraph *g)
{
    for (VertexIdx i=0; i < g->nVertices; ++i)
        fprintf(f, "%ld: %ld\n",i,info->perVertex[i]);
    printf("---------------------\n");
    for (VertexIdx i=0; i < g->nVertices; ++i)
        for (EdgeIdx j=g->offsets[i]; j < g->offsets[i+1]; ++j)
            fprintf(f, "(%ld, %ld): %ld\n",i,g->nbors[j],info->perEdge[j]);
}


#endif



