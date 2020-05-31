import numpy as np
import sys
import itertools
import pandas as pd
from keras.models import Sequential
from keras import layers
from matplotlib import pyplot
from sklearn.datasets import make_circles
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import roc_auc_score

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
    hist = model.fit(X, y, epochs=50, batch_size=10, verbose=0)
    predictions = model.predict_classes(X, verbose=0)
    # accuracy: (tp + tn) / (p + n)
    accuracy = accuracy_score(y, predictions)
    print('Accuracy ((TP + TN) / (P + N)): %f' % accuracy)
    # precision tp / (tp + fp)
    precision = precision_score(y, predictions)
    print('Precision (TP / (TP + FP)): %f' % precision)
    # recall: tp / (tp + fn)
    recall = recall_score(y, predictions)
    print('Recall (TP / (TP + FN)): %f' % recall)
    # f1: 2 tp / (2 tp + fp + fn)
    f1 = f1_score(y, predictions)
    print('F1 Score (2TP / (2TP + FP + FN)): %f' % f1)
