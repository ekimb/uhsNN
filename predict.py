import argparse
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import load_model
import numpy as np
import pandas as pd
import os
import itertools
import sys
from timeit import default_timer as timer
from multiprocessing import Pool
def buildDict(k):
    bases = ['A', 'C', 'G', 'T']
    basedict={"A": 0, "C": 1, "G": 2, "T": 3}
    intdict = {}
    kmerArray = [''.join(p) for p in itertools.product(bases, repeat=k)]
    for kmer in kmerArray:
        totval = 0
        for pos,char in enumerate(kmer):
            val = basedict[char] * pow(4, (k - pos - 1))
            totval += val
        intdict[totval] = kmer
    return intdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Makes predictions for provided k-mer length (k) and sequence length (L) with trained model (.h5) as input.')
    parser.add_argument('-k', metavar='k', type=int, help='K-mer size (k) for the UHS')
    parser.add_argument('-L', metavar= 'L', type=int, help='Sequence size (L) for the UHS')
    parser.add_argument('-m', metavar='m', help='Complete model (.h5) file')
    parser.add_argument('-t', metavar='t', help='Model threshold')
    parser.add_argument('-n', metavar= 'n', help='Number of threads')


    args = parser.parse_args()
    k = args.k
    L = args.L
    m = args.m
    decycling=pd.read_csv("int/decyc"+str(k)+".int")
    maxmaxk = 15
    intdict = buildDict(k)
    model = load_model(args.m, compile=False)
    print("Loaded model from disk")
    print(L)
    lst = np.array(list(itertools.product([0, 1], repeat=2*k)))
    padlst = np.pad(lst, ((0, 0),(0,maxmaxk*2-lst.shape[1])), 'constant', constant_values=((2, 2),(2,2)))
    lstL = np.append(np.full(lst.shape[0], L).reshape(lst.shape[0],1), padlst, axis=1)
    lstkL = np.append(np.full(lst.shape[0], k).reshape(lstL.shape[0],1), lstL, axis=1)
    start = timer()
    labels=model.predict(lstkL, batch_size=8192, workers=args.n, use_multiprocessing=True)
    keras.backend.learning_phase(0)
    end = timer()
    labels[decycling]=2
    print(labels.shape)
    outF = open("preds.txt", "w")
    for i in range(len(labels)):
        outF.write(intdict[i] + '\t' + str(labels[i][0]))
        outF.write("\n")
    outF.close()
    predTime = end - start
    print("Predictions done, took " + str(predTime) + " seconds.")
    print("Starting DOCKS...")
    start = timer()
    command = "./predict preds.txt " + str(k) + " " + str(L) + " decyc" + str(k) + ".txt " +  str(k) + str(L) + "_predicted.txt " + str(args.t) + " " + str(args.n)
    print(command)
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname + "/build/")
    stream = os.popen(command)
    output = stream.read()
    print(output)
    end = timer()
    DOCKStime = end - start
    print("PDOCKS done, took " + str(DOCKStime) + " seconds.")
    print("Total is " + str(predTime+DOCKStime) + " seconds.")