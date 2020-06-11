#include <iostream>
#include <sys/resource.h> 
#include <sys/time.h> 
#include <unistd.h> 
#include <stdio.h> 
#include <time.h>
#include<bits/stdc++.h>
#include <sys/stat.h> 
#include <sys/types.h> 
#include "src/decycling.h"
#include "src/docks.h"

int main(int argc, char* argv[]) {
    // Initialize model.
    string modelPath = argv[1];
    // Create a 1D Tensor on length 10 for input data.
    struct timespec start, finish;
    double elapsed;
    size_t pos;
    int k = stoi(argv[2], &pos);
    int L = stoi(argv[3], &pos);
    const char *decyclingPath = argv[4];
    const char *hittingPath = argv[5];
    double threshold = stod(argv[6]);
    int threads = stoi(argv[6]);
    const double PI = 3.14159;
    string directory;
    ofstream decyclingStream(decyclingPath);
    DOCKS docks = DOCKS(k);
    decycling newDecycling;
    vector<unsigned int> decyclingSet = newDecycling.computeDecyclingSet(k);
    for (int i = 0; i < decyclingSet.size(); i++) {
        string label = docks.getLabel(decyclingSet[i]);
        docks.removeEdge(decyclingSet[i]);
        decyclingStream << label << "\n";
    }
    int decyclingSize = decyclingSet.size();
    cout << "Decycling set size: " << decyclingSize << endl;
    decyclingStream.close();
    //clock_gettime(CLOCK_MONOTONIC, &start);
    int hittingSize = docks.HittingML(L, hittingPath, modelPath, threshold, threads);
    //clock_gettime(CLOCK_MONOTONIC, &finish);
    //elapsed = (finish.tv_sec - start.tv_sec);
    //elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    cout << elapsed << " seconds." << endl; 
    cout << "Hitting set size: " << hittingSize << endl << endl;
    return 0;
}
