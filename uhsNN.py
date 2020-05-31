import numpy as np
import sys
import itertools
import pandas as pd
from keras.models import Sequential
from keras import layers

from sklearn.model_selection import train_test_split


def oneHotEncode(kmer):
    mapping = dict(zip("ACGT", range(4)))    
    seq2 = [mapping[i] for i in kmer]
    return np.eye(4)[seq2]

def encodeReverse(kmerList):
    kmer = ""
    for i in range(len(kmerList)):
        if kmerList[i] == 1:
            if i % 4 == 0:
                kmer = kmer + 'A'
            elif i % 4 == 1:
                kmer = kmer + 'C'
            elif i % 4 == 2:
                kmer = kmer + 'G'
            elif i % 4 == 3:
                kmer = kmer + 'T'
    return kmer
if __name__ == "__main__":
    bases = ['A', 'C', 'G', 'T']
    UHSpath = sys.argv[1]
    k = int(sys.argv[2])
    kmerArray = [''.join(p) for p in itertools.product(bases, repeat=k)]
    f = open(UHSpath, "r")
    UHSarray = []
    trainInput = []
    trainOutput = []
    for line in f:
        UHSkmer = line.rstrip()
        UHSarray.append(UHSkmer)
    for kmer in kmerArray:
        trainInput.append(np.array(oneHotEncode(kmer).flatten()))
        if kmer in UHSarray:
            trainOutput.append(1)
        else:
            trainOutput.append(0)
    model = Sequential()
    model.add(layers.Dense(4, input_dim=4*k, activation='relu'))
    model.add(layers.Dense(4, activation='relu'))
    model.add(layers.Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    model.summary()
    X = np.array(trainInput)
    y = np.array(trainOutput)
    model.fit(X, y, epochs=50, batch_size=10)
    _, accuracy = model.evaluate(X, y)
    print('Accuracy: %.2f' % (accuracy*100)) 
