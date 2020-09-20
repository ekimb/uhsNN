import argparse
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import load_model
import numpy as np
import os
import itertools
import sys
from train import enumerateKmers, oneHotEncode, constructUHSarray
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
    parser.add_argument('-L', metavar= 'L', help='Sequence size (L) for the UHS')
    parser.add_argument('-m', metavar='m', help='Complete model (.h5) file')
    parser.add_argument('-t', metavar='t', help='Model threshold')
    parser.add_argument('-n', metavar= 'n', help='Number of threads')


    start = timer()
    args = parser.parse_args()
    k = args.k
    L = args.L
    m = args.m
    maxmaxk = 15
    intdict = buildDict(k)
    model = load_model(args.m, compile=False)
    print("Loaded model from disk")
    print(L)
    total_labels=np.empty(0,dtype=np.int8)
    total_lst=np.empty([0,2*maxmaxk+2],dtype=np.int8)
    lst = np.array(list(itertools.product([0, 1], repeat=2*k)))
    padlst = np.pad(lst, ((0, 0),(0,maxmaxk*2-lst.shape[1])), 'constant', constant_values=((2, 2),(2,2)))
    lstL = np.append(np.full(lst.shape[0], L).reshape(lst.shape[0],1), padlst, axis=1)
    lstkL = np.append(np.full(lst.shape[0], k).reshape(lstL.shape[0],1), lstL, axis=1)
    labels=np.zeros(np.power(4,k), dtype=np.int8)
    print(labels.shape)
    combined = np.append(lstkL, labels.reshape(labels.shape[0],1), axis=1)
    total_lst=np.append(total_lst, combined[:,0:2*maxmaxk+2],axis=0)
    prediction = model.predict(total_lst)
    outF = open("preds.txt", "w")
    for i in range(len(prediction)):
        outF.write(intdict[i] + '\t' + str(prediction[i][0]))
        outF.write("\n")
    outF.close()
    end = timer()
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
    print("DOCKS done, took " + str(DOCKStime) + " seconds.")
    print("Total is " + str(predTime+DOCKStime) + " seconds.")