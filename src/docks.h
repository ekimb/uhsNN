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
using keras2cpp::Model;
using keras2cpp::Tensor;
using namespace std;
using byte8 = uint8_t;
class DOCKS {
    public:
        byte8* finished;
        byte8* used;
        double* hittingNumArray;
        float** D;
        float** F;
        int ALPHABET_SIZE;
        double edgeCount;
        double edgeNum;
        int k;
        int curr;
        int l;
        int total;
        int vertexCount; 
        int vertexExp;
        int vertexExp2;
        int vertexExp3;
        int vertexExpMask;
        int vertexExp_1;
        byte8* edgeArray;
        int* topoSort;
        static map<char, int> alphabetMap;
        string ALPHABET;
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
        edgeArray = new byte8[(int)edgeNum];
        generateGraph(k);
        map<char, int> alphabetMap;
        for (int i = 0; i < ALPHABET_SIZE; i++) alphabetMap.insert(pair<char,int>(ALPHABET[i], i));
    }

    void generateGraph(int k) {  
    /**
    Generates a complete de Bruijn graph of order k.
    @param k: Desired k-mer length (order of complete graph).
    */
        for (int i = 0; i < edgeNum; i++) edgeArray[i] = 1;
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
                if ((depth[vertexIndex] > maxVertDepth) && (edgeArray[edgeIndex] == 1)) maxVertDepth = depth[vertexIndex];
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
        if (edgeArray[i] == 1) edgeCount--;
        edgeArray[i] = 0;
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
            if (edgeArray[index] == 1) adjVertex[count++] = index / ALPHABET_SIZE;
        }
        vector<int> rc(count);
        for (int i = 0; i < count; i++) {
            rc[i] = adjVertex[i];
        }
        return rc;
    }

    int Hitting(int L, const char *hittingPath) {
    /**
    Performs hitting set calculations with parallelization
    and without randomization, counting L-k+1-long paths.
    @param L: Sequence length, hittingFile: Output file destination.
    @return hittingCount: Size of hitting set.
    */
        vertexExp = pow(ALPHABET_SIZE, k-1);
        int imaxHittingNum = -1;
        ofstream hittingStream;
        int hittingCount = 0;
        l = L-k+1;
        hittingNumArray = new double[(int)edgeNum];
        used = new byte8[vertexExp];
        finished = new byte8[vertexExp];
        topoSort = new int[vertexExp];
        D = new float*[l + 1];
        float* Dpool = new float[(l+1)* vertexExp];
        for(int i = 0; i < l+1; i++, Dpool += vertexExp) D[i] = Dpool;
        hittingStream.open(hittingPath); 
        F = new float*[l + 1];
        float* Fpool = new float[(l+1)* vertexExp];
        for(int i = 0; i < l+1; i++, Fpool += vertexExp) F[i] = Fpool;
        while (calculatePaths(l)) {
            int imaxHittingNum = calculateHittingNumber(l);
            if (imaxHittingNum < 0) break;
            removeEdge(imaxHittingNum);
            string label = getLabel(imaxHittingNum);
            hittingStream << label << "\n";
            hittingCount++;
        }
        hittingStream.close();
        topologicalSort();
        cout << "Length of longest remaining path: " <<  maxLength() << "\n";
        return hittingCount;
    }

    int calculateHittingNumber(int L) {
/**
Calculates hitting number of all edges, counting paths of length L-k+1, in parallel.
@param L: Sequence length.
@return imaxHittingNum: Index of vertex with maximum hitting number.
*/  
        double maxHittingNum = 0;
        int imaxHittingNum = -1;
        for (int i = 0; i < edgeNum; i++) {
            double hittingNum = 0;
            for (int j = (1 - edgeArray[i]) * L; j < L; j++) {
                hittingNum = hittingNum + F[j][i % vertexExp]/1.4e-45 * D[(L-j-1)][i / ALPHABET_SIZE]/1.4e-45;
                hittingNumArray[i] = hittingNum;
            }
            if (hittingNum > maxHittingNum) {maxHittingNum = hittingNum; imaxHittingNum = i;}
        }
        return imaxHittingNum;
    }
    Tensor makePrediction (int index, Model &model) {
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
        return out;
    }
    int HittingML(int L, const char *hittingPath, string modelPath, double threshold) {
    /**
    Performs hitting set calculations with parallelization
    and without randomization, counting L-k+1-long paths.
    @param L: Sequence length, hittingFile: Output file destination.
    @return hittingCount: Size of hitting set.
    */
        auto model = Model::load(modelPath);
        vertexExp = pow(ALPHABET_SIZE, k-1);
        int imaxHittingNum = -1;
        ofstream hittingStream;
        int hittingCount = 0;
        l = L-k+1;
        hittingNumArray = new double[(int)edgeNum];
        used = new byte8[vertexExp];
        finished = new byte8[vertexExp];
        topoSort = new int[vertexExp];
        D = new float*[l + 1];
        float* Dpool = new float[(l+1)* vertexExp];
        for(int i = 0; i < l+1; i++, Dpool += vertexExp) D[i] = Dpool;
        hittingStream.open(hittingPath); 
        F = new float*[l + 1];
        float* Fpool = new float[(l+1)* vertexExp];
        for(int i = 0; i < l+1; i++, Fpool += vertexExp) F[i] = Fpool;
        std::cout << "here" << std::endl;
        for (int i = 0; i < edgeNum; i++) {
            if (edgeArray[i] == 1) {
                Tensor tout = makePrediction(i, model);
                double res;
                for (auto&& it : tout.data_) {
                    res = static_cast<double>(it);
                }
                if (res >= threshold) {
                    removeEdge(i);
                    string label = getLabel(i);
                    hittingStream << label << "\n";
                    hittingCount++;
                }
            }
        }
        topologicalSort();
        cout << "Length of longest remaining path after model prediction: " <<  maxLength() << "\n";
        while (calculatePaths(l) == 1) {
            int imaxHittingNum = calculateHittingNumber(l);
            if (imaxHittingNum < 0) break;
            removeEdge(imaxHittingNum);
            string label = getLabel(imaxHittingNum);
            hittingStream << label << "\n";
            hittingCount++;
        }
        hittingStream.close();
        topologicalSort();
        cout << "Length of longest remaining path: " <<  maxLength() << "\n";
        return hittingCount;
    }

    int calculatePaths(int L) {
    /**
    Calculates number of L-k+1 long paths for all vertices.
    @param L: Sequence length.
    @return 1: True if path calculation completes.
    */
        vertexExp2 = vertexExp * 2;
        vertexExp3 = vertexExp * 3;
        vertexExpMask = vertexExp - 1;
        vertexExp_1 = pow(ALPHABET_SIZE, k-2);
        for (int i = 0; i < vertexExp; i++) {D[0][i] = 1.4e-45; F[0][i] = 1.4e-45;}
        for (int j = 1; j <= L; j++) {
            for (int i = 0; i < vertexExp; i++) {
                int index = (i * 4);
                F[j][i] = edgeArray[index]*F[j-1][index & vertexExpMask] + edgeArray[index + 1]*F[j-1][(index + 1) & vertexExpMask] + edgeArray[index + 2]*F[j-1][(index + 2) & vertexExpMask] + edgeArray[index + 3]*F[j-1][(index + 3) & vertexExpMask];
                D[j][i] = edgeArray[i]*D[j-1][(i >> 2)] + edgeArray[i + vertexExp]*D[j-1][((i + vertexExp) >> 2)] + edgeArray[i + vertexExp2]*D[j-1][((i + vertexExp2) >> 2)] + edgeArray[i + vertexExp3]*D[j-1][((i + vertexExp3) >> 2)];
            }
        }
    return 1;
    }
};
#endif