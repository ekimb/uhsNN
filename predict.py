import argparse
from keras.models import Sequential
from keras.layers import Dense
from keras.models import model_from_json
import numpy as np
import os
import itertools
import sys
from train import enumerateKmers, oneHotEncode, constructUHSarray
from timeit import default_timer as timer
from multiprocessing import Pool

def constructDecycArray(path):
    f = open(path, "r")
    decycArray = []
    for line in f:
        if line[0] == '>decycling':
            decycKmerSwitch = 1
            addKmerSwitch = 0
        elif line[0] == '>':
            break
        else:
            kmer = line.rstrip()
            decycArray.append(kmer)
    return decycArray

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trains a model with additional k-mer sets (.txt files) as input.')
    parser.add_argument('-k', metavar='k', type=int, help='K-mer size (k) for the UHS')
    parser.add_argument('-L', metavar= 'L', help='Sequence size (L) for the UHS')
    parser.add_argument('-u', metavar='u', help='UHS (.uhs) file')
    parser.add_argument('-m', metavar='m', help='Model (.model) file')
    parser.add_argument('-t', metavar='t', help='Model threshold')
    parser.add_argument('-n', metavar= 'n', help='Number of threads')

    start = timer()
    args = parser.parse_args()
    k = args.k
    L = args.L
    UHSfile = args.u
    jsonFile = open('model8.json', 'r')
    loadedModelJson = jsonFile.read()
    jsonFile.close()
    loadedModel = model_from_json(loadedModelJson)
    # load weights into new model
    loadedModel.load_weights("model8.h5")
    print("Loaded model from disk")
    bases = ['A', 'C', 'G', 'T']
    kmerArray = [''.join(p) for p in itertools.product(bases, repeat=k)]
    decycArray = constructDecycArray(UHSfile)
    predictInput = []
    kmerArray = [x for x in kmerArray if x not in decycArray]
    print('Preparing k-mers for k = %d and L = %s' % (k, L))
    for kmer in kmerArray:
        oneHot = oneHotEncode(kmer).flatten()
        encodedKmer = np.append(oneHot, int(L))
        predictInput.append(np.array(encodedKmer))
    prediction = loadedModel.predict(np.array(predictInput))
    outF = open("preds.txt", "w")
    for i in range(len(prediction)):
        outF.write(kmerArray[i] + '\t' + str(prediction[i][0]))
        outF.write("\n")
    outF.close()
    end = timer()
    predTime = end - start
    print("Predictions done, took " + str(predTime) + " seconds.")
    print("Starting DOCKS...")
    start = timer()
    command = "./predict preds.txt " + str(k) + " " + str(L) + " decyc8.txt " +  "8105.txt " + str(args.t) + " " + str(args.n)
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