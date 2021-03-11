#ifndef DOCKS_H
#define DOCKS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cstdlib>
#include <iomanip>
#include <cstdint>
#include <omp.h>
#include "model.h"
#include "tensor.h"
#include "utils.h"
#include <stdlib.h>
using keras2cpp::Model;
using keras2cpp::Tensor;
using keras2cpp::Stream;
using namespace std;
using byte8 = uint8_t;
struct record
{
   int index;
   std::string kmer;
   float pred;
};
class DOCKS {
    public:
        byte8* finished;
        byte8* used;
        double* hittingNumArray;
        float** D;
        float* Fcurr;
        float delta;
        float epsilon;
        float* Fprev;
        int ALPHABET_SIZE;
        double edgeCount;
        double edgeNum;
        int k;
        int curr;
        int l;
        int h;
        int total;
        int exit;
        int vertexCount; 
        int vertexExp;
        int vertexExp2;
        int vertexExp3;
        int vertexExpMask;
        int vertexExp_1;
        vector<bool> edgeArray;
        byte8* stageArray;
        byte8* pick;
        int* topoSort;
        static map<char, int> alphabetMap;
        string ALPHABET;
        vector<int> stageVertices;
    DOCKS (int argK) {
    /**
    Definition of a graph object. Generates a graph of order k, creates an empty
    edge index array, calculates number of edges, builds a character-index map.
    @param argK: Argument passed as k-mer length.
    */
        
        ALPHABET = "ACGT";
        ALPHABET_SIZE = 4;
        k = argK;
        edgeNum = pow(ALPHABET_SIZE, k);
        cout << "Constructing edges..." << endl;
        vector<bool> edgeArray(edgeNum, true);
        cout << "Constructed edges."  << endl;
        generateGraph(k);
        map<char, int> alphabetMap;
        for (int i = 0; i < ALPHABET_SIZE; i++) alphabetMap.insert(pair<char,int>(ALPHABET[i], i));
    }

    void generateGraph(int k) {  
    /**
    Generates a complete de Bruijn graph of order k.
    @param k: Desired k-mer length (order of complete graph).
    */
        for (int i = 0; i < edgeNum; i++) edgeArray.push_back(true);
        edgeCount = edgeNum;
        vertexCount = edgeNum / ALPHABET_SIZE; 
    }

    char getChar(int i) {
    /**
    Gets alphabet character from index.
    @param i: Index of character.
    @return The character in the alphabet.
    */
        return ALPHABET[i];
    }

    string getLabel(int i) {
    /**
    Gets label of the input edge index.
    @param i: Index of edge.
    @return The label of the edge.
    */
        string finalString = "";
        for (int j = 0; j < k; j++) {
            finalString = getChar((i % ALPHABET_SIZE)) + finalString;
            i = i / ALPHABET_SIZE;
        }
        return finalString;
    }
    int getIndex(string label) {
    /**
    Gets index of the input edge label.
    @param i: label of edge.
    @return The index of the edge.
    */
        map<char, int> alphabetMap;
        for (int i = 0; i < ALPHABET_SIZE; i++) alphabetMap.insert(pair<char,int>(ALPHABET[i], i));
        int index = 0;
        for (int j = 0; j < k; j++) {
            //cout << alphabetMap[label[j]] << endl;
            index += alphabetMap[label[j]] * pow(4, k-j-1);
            //cout << index << endl;
        }
        return index;
    }
    int maxLength() {
    /**
    Calculates the length of the maximum length path in the graph.
    @return maxDepth: Maximum length.
    */
        vector<int> depth(vertexExp);
        int maxDepth = -1;
        for (int i = 0; i < vertexExp; i++) {
            int maxVertDepth = -1;
            for (int j = 0; j < ALPHABET_SIZE; j++) {
                int edgeIndex = topoSort[i] + j * vertexExp;
                int vertexIndex = edgeIndex / ALPHABET_SIZE;
                if ((depth[vertexIndex] > maxVertDepth) && (edgeArray[edgeIndex] == true)) maxVertDepth = depth[vertexIndex];
            }
            depth[topoSort[i]] = maxVertDepth + 1;
            if (depth[topoSort[i]] > maxDepth) {maxDepth = depth[topoSort[i]];}
        }
        return maxDepth;
    }

    void removeEdge(int i) {
    /**
    Removes an edge from the graph.
    @param i: Index of edge.
    */
        if (edgeArray.at(i) == true) edgeCount--;
        edgeArray.at(i) = false;
    }

    void topologicalSort() {
    /**
    Traverses the graph in topological order.
    */
        for (int i = 0; i < vertexExp; i++) {used[i] = false; finished[i] = false;}
        int index = 0;
        for (int i = 0; i < vertexExp; i++) {
            if (used[i] == false) {
                index = depthFirstSearch(index, i);
                if (index == -1) {topoSort = NULL; return;}
            }
        }
       // int rc[vertexExp];
       // for (int i = 0; i < vertexExp; i++) rc[i] = topoSort[vertexExp-i-1];
    }
    int depthFirstSearch(int index, int u) {
    /**
    Depth-first search of a given index of an edge.
    @param index: Depth of recursion, u: Index of edge.
    @return -1: The search cycles, index+1: Current depth.
    */
        used[u] = true;
        bool cycle = false;
        for (int v : getAdjacent(u)) {
            if (used[v] == true && finished[v] == false) cycle = true;
            if (used[v] == false) {
                index = depthFirstSearch(index, v);
                cycle = cycle || (index == -1);
            }
        }
        finished[u] = true;
        topoSort[index] = u;
        if (cycle) return -1;
        else return index + 1;
    }
    vector<int> getAdjacent(int v) {
    /**
    Get adjacent vertices to a given index of a vertex.
    @param v: Index of vertex.
    @return rc: Array of adjacent vertices.
    */
        int count = 0;
        int adjVertex[ALPHABET_SIZE];
        for (int i = 0; i < ALPHABET_SIZE; i++) {
            int index = v + i * vertexExp;
            if (edgeArray[index] == true) adjVertex[count++] = index / ALPHABET_SIZE;
        }
        vector<int> rc(count);
        for (int i = 0; i < count; i++) {
            rc[i] = adjVertex[i];
        }
        return rc;
    }


    int calculatePaths(int L, int threads) {
    /**
    Calculates number of L-k+1 long paths for all vertices.
    @param L: Sequence length.
    @return 1: True if path calculation completes.
    */
        omp_set_dynamic(0);
        curr = 1;
        vertexExp2 = vertexExp * 2;
        vertexExp3 = vertexExp * 3;
        vertexExpMask = vertexExp - 1;
        vertexExp_1 = pow(ALPHABET_SIZE, k-2);
        int CHUNK_SIZE = vertexExp / threads;
        #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(threads)
        for (int i = 0; i < vertexExp; i++) {D[0][i] = 1.4e-45;Fprev[i] = 1.4e-45;}
        for (int j = 1; j <= L; j++) {
            #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(threads)
            for (int i = 0; i < vertexExp; i++) {
                D[j][i] = edgeArray[i]*D[j-1][(i >> 2)] + edgeArray[i + vertexExp]*D[j-1][((i + vertexExp) >> 2)] + edgeArray[i + vertexExp2]*D[j-1][((i + vertexExp2) >> 2)] + edgeArray[i + vertexExp3]*D[j-1][((i + vertexExp3) >> 2)];
            }
        }
        int CHUNK_SIZE_EDGE = (int)edgeNum / threads;
        #pragma omp parallel for schedule(static, CHUNK_SIZE_EDGE) num_threads(threads)
        for (int i = 0; i < (int)edgeNum; i++) hittingNumArray[i] = 0;
        while (curr <= L) {
            #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(threads)
            for (int i = 0; i < vertexExp; i++) {
                int index = (i * 4);
                Fcurr[i] = edgeArray[index]*Fprev[index & vertexExpMask] + edgeArray[index + 1]*Fprev[(index + 1) & vertexExpMask] + edgeArray[index + 2]*Fprev[(index + 2) & vertexExpMask] + edgeArray[index + 3]*Fprev[(index + 3) & vertexExpMask];
            }
            #pragma omp parallel for schedule(static, CHUNK_SIZE_EDGE) num_threads(threads)
            for (int i = 0; i < (int)edgeNum; i++) {
                hittingNumArray[i] += (Fprev[i % vertexExp]/1.4e-45) * (D[(L-curr)][i / ALPHABET_SIZE]/1.4e-45);
                if (edgeArray[i] == false) hittingNumArray[i] = 0;
            }
            #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(threads)
            for (int i = 0; i < vertexExp; i++) Fprev[i] = Fcurr[i];
            #pragma omp single
            curr++;
        }
        return 1;
    }
    Tensor makePrediction (int index, Model &model, int threads) {
        Tensor t{4*k+1};
        string kmer = getLabel(index);
        for (int i = 0; i < k; i++) {
            if (kmer[i] == 'A') {
                t.data_[4*i] = 1;
                t.data_[4*i+1] = 0;
                t.data_[4*i+2] = 0;
                t.data_[4*i+3] = 0;
            }
            else if (kmer[i] == 'C') {
                t.data_[4*i] =0;
                t.data_[4*i+1] = 1;
                t.data_[4*i+2] = 0;
                t.data_[4*i+3] = 0;
            }
            else if (kmer[i] == 'G') {
                t.data_[4*i] = 0;
                t.data_[4*i+1] = 0;
                t.data_[4*i+2] = 1;
                t.data_[4*i+3] = 0;
            }
            else if (kmer[i] == 'T') {
                t.data_[4*i] = 0;
                t.data_[4*i+1] = 0;
                t.data_[4*i+2] = 0;
                t.data_[4*i+3] = 1;
            }
        }
        t.data_[4*k] = (l+k-1);
        Tensor out = model(t);
        t.data_.clear();
        t.data_.shrink_to_fit();
        return out;
    }
    int HittingRandomParallel(int L, const char *hittingPath, float threshold, int threads, vector<record> &v) {
    /**
    Performs hitting set calculations with parallelization
    and with randomization, counting L-k+1-long paths.
    @param L: Sequence length, hittingFile: Output file destination.
    @return hittingCount: Size of hitting set.
    */
        omp_set_dynamic(0);
        srand (1);
        int hits = 0;
        vertexExp = pow(ALPHABET_SIZE, k-1);
        ofstream hittingStream(hittingPath);
        int hittingCount = 0;
        l = L-k+1;
        int i, j;
        delta = 1/(double)l;
        epsilon = (1-8*(delta))/4;
        double alpha = 1 - 4*delta -2*epsilon;
        cout << "Alpha: " << 1/alpha << endl << "Delta: " << delta << endl << "Epsilon: " << epsilon << endl;
        hittingNumArray = new double[(int)edgeNum];
        stageArray = new byte8[(int)edgeNum];
        used = new byte8[vertexExp];
        finished = new byte8[vertexExp];
        pick = new byte8[(int)edgeNum];
        topoSort = new int[vertexExp];
        D = new float*[l + 1];
        Fcurr = new float[vertexExp];
        Fprev = new float[vertexExp];
        #pragma omp parallel for schedule(static) num_threads(threads)
        for(int i = 0; i < l+1; i++) {D[i] = new float[vertexExp];}
        topologicalSort();
        cout << "Length of longest remaining path after model prediction: " <<  maxLength() << "\n";
        calculatePaths(l, threads);
        int imaxHittingNum = calculateHittingNumberParallel(l, false, threads);
        cout << "Max hitting number: " << hittingNumArray[imaxHittingNum] << endl;
        h = findLog((1.0+epsilon), hittingNumArray[imaxHittingNum]);
        double prob = delta/(double)l;
        while (h > 0) {
            cout << "H is " << h << endl;
            total = 0;
            int hittingCountStage = 0;
            double pathCountStage = 0;
            calculatePaths(l, threads);
            imaxHittingNum = calculateHittingNumberParallel(l, true, threads);
            cout << "Calculated hitting numbers." << endl;
            if (exit == -1) break;
            stageVertices = pushBackVector(threads);
            cout << "Stage size: " << stageVertices.size() << endl;
            int CHUNK_SIZE = stageVertices.size() / threads;
            if (stageVertices.size() <= threads) {
                CHUNK_SIZE = 1;
            }
            #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(threads)
            for (int it = 0; it < stageVertices.size(); it++) {
                i = stageVertices[it];
                if ((!pick[i]) && (hittingNumArray[i] > (pow(delta, 3) * total))) {
                    stageArray[i] = 0;
                    pick[i] = true;
                    #pragma omp critical
                    {
                        hittingCountStage++;
                        pathCountStage += hittingNumArray[i];
                    }
                }
            }
            #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(threads)
            for (int it = 0; it < stageVertices.size(); it++) {
                for (int jt = 0; jt < stageVertices.size(); jt++) {
                    i = stageVertices[it];
                    if (!pick[i]) {
                        if (((double) rand() / (RAND_MAX)) <= prob) {
                            stageArray[i] = 0;
                            pick[i] = true;
                            #pragma omp critical
                            {
                                hittingCountStage += 1;
                                pathCountStage += hittingNumArray[i];
                            }
                        }
                        j = stageVertices[jt];
                        if (!pick[j]) {
                            if (((double) rand() / (RAND_MAX)) <= prob) {
                                stageArray[j] = 0;
                                pick[j] = true;
                                #pragma omp critical
                                {
                                    hittingCountStage += 1;
                                    pathCountStage += hittingNumArray[i];
                                }

                            }
                            else pick[i] = false;
                        }
                    }
                }
            }
            hittingCount += hittingCountStage;
            if (pathCountStage >= hittingCountStage * pow((1.0 + epsilon), h) * (1 - 4*delta - 2*epsilon)) {
                #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(threads)
                for (int it = 0; it < stageVertices.size(); it++) {
                    i = stageVertices[it];
                    if (pick[i] == true) {
                        removeEdge(i);
                        string label = getLabel(i);
                        #pragma omp critical 
                        {
                            hittingStream << label << "\n";
                            hits++;
                        }
                    }
                }
                h--;
            }
            else {
                hittingCount -= hittingCountStage;
            }
        }
        hittingStream.close();
        topologicalSort();
        cout << "Length of longest remaining path: " <<  maxLength() << "\n";
        return hits;
        
    }
    int findLog(double base, double x) {
    /**
    Finds the logarithm of a given number with respect to a given base.
    @param base: Base of logartihm, x: Input number.
    @return (int)(log(x) / log(base)): Integer logarithm of the number and the given base.
    */
        return (int)(log(x) / log(base));
    }
    vector<int> pushBackVector(int threads) {
        vector<int> stageVertices;
        int CHUNK_SIZE = (int)edgeNum / threads;
        for(int i = 0; i < (int)edgeNum; i++) {
            if (stageArray[i] == 1) stageVertices.push_back(i);
        }
        return stageVertices;
    }
    int calculateHittingNumberParallel(int L, bool random, int threads) {
/**
Calculates hitting number of all edges, counting paths of length L-k+1, in parallel.
@param L: Sequence length.
@return imaxHittingNum: Index of vertex with maximum hitting number.
*/  
        omp_set_dynamic(0);
        double maxHittingNum = 0;
        int imaxHittingNum = 0;
        int count = 0;
        exit = -1;
        int CHUNK_SIZE = (int)edgeNum / threads;
        #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(threads)
        for (int i = 0; i < (int)edgeNum; i++) {
            if (random == true) {
                if (((hittingNumArray[i]) >= pow((1.0+epsilon), h-1)) && ((hittingNumArray[i]) <= pow((1.0+epsilon), h))) {
                    stageArray[i] = 1;
                    pick[i] = false;
                    #pragma omp critical
                    {
                        total += hittingNumArray[i] * stageArray[i];
                        count++;
                    }
                }
                else {
                    stageArray[i] = 0;
                    pick[i] = false;
                }   
            }
        }
        #pragma omp parallel for schedule(static, CHUNK_SIZE) num_threads(threads)
        for (int i = 0; i < (int)edgeNum; i++) {
            if ((hittingNumArray[i])*edgeArray[i] > maxHittingNum) {
                maxHittingNum = hittingNumArray[i]; imaxHittingNum = i;
                exit = 0;
            }
        }
        return imaxHittingNum;
    }
};
#endif