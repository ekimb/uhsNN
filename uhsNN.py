import numpy as np
import sys
import itertools
import pandas as pd
import tensorflow as tf
import argparse
from keras.models import Sequential
from keras import layers
from matplotlib import pyplot
from sklearn.datasets import make_circles
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

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

def constructUHSarray(path):
    f = open(path, "r")
    decycArray = []
    addArray = {}
    decycKmerSwitch = 0
    addKmerSwitch = 0
    for line in f:
        if line[0] == '>' and line.strip('\n').strip('>') == 'decycling':
            decycKmerSwitch = 1
            addKmerSwitch = 0
        elif line[0] == '>':
            decycKmerSwitch = 0
            addKmerSwitch = 1
            seqLength = int(line.strip('\n').strip('>'))
            addArray[seqLength] = []
        else:
            if decycKmerSwitch == 1 and addKmerSwitch == 0:
                kmer = line.rstrip()
                decycArray.append(kmer)
            elif decycKmerSwitch == 0 and addKmerSwitch == 1:
                kmer = line.rstrip()
                addArray[seqLength].append(kmer)
    return decycArray, addArray
def enumerateKmers(k):
    bases = ['A', 'C', 'G', 'T']
    kmerArray = [''.join(p) for p in itertools.product(bases, repeat=k)]
    return kmerArray
def decyclingPredict(decycPath, k):
    trainInput = []
    trainOutput = []
    UHSarray = constructUHSarray(decycPath)
    kmerArray = enumerateKmers(k)
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
    model.fit(X, y, epochs=20, batch_size=10)
    predictions = model.predict_classes(X, verbose=0)
    accuracy = accuracy_score(y, predictions)
    print('Accuracy ((TP + TN) / (P + N)): %f' % accuracy)
    precision = precision_score(y, predictions)
    print('Precision (TP / (TP + FP)): %f' % precision)
    recall = recall_score(y, predictions)
    print('Recall (TP / (TP + FN)): %f' % recall)
    f1 = f1_score(y, predictions)
    print('F1 Score (2TP / (2TP + FP + FN)): %f' % f1)

def additionalPredict(k, UHSfile, epochs, batchSize):
    trainInput = []
    trainOutput = []
    decycArray, addArray = constructUHSarray(UHSfile)
    kmerArray = enumerateKmers(k)
    for kmer in kmerArray:
        if kmer in decycArray:
            kmerArray.remove(kmer)
    for L in addArray.keys():
        print('Reading k-mers for k = %d and L = %s' % (k, L))
        for kmer in kmerArray:
            encodedKmer = np.append(oneHotEncode(kmer).flatten(), int(L))
            trainInput.append(np.array(encodedKmer))
            if kmer in addArray[L]:
                trainOutput.append(1)
            else:
                trainOutput.append(0)
    model = Sequential()
    model.add(layers.Dense(64, input_dim=4*k+1, activation='relu'))
    model.add(layers.Dense(8, activation='relu'))
    model.add(layers.Dense(8, activation='relu'))
    model.add(layers.Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=[tf.keras.metrics.PrecisionAtRecall(recall=0.99)])
    model.summary()
    X = np.array(trainInput)
    y = np.array(trainOutput)
    addModel = model.fit(X, y, epochs=epochs, batch_size=batchSize)
    with open('./modelArch.json', 'w') as fout:
        fout.write(model.to_json())
    model.save_weights('./modelWeights.h5', overwrite=True)
    predictions = model.predict_proba(X)
    auc = roc_auc_score(y, predictions)
    print("AUC: ", auc)
    fpr, tpr, thresholds = roc_curve(y, predictions)

    # plot no skill
    pyplot.plot([0, 1], [0, 1], linestyle='--')
    # plot the roc curve for the model
    pyplot.plot(fpr, tpr)
    # show the plot
    pyplot.show()


    '''for i in range(len(predictions)):
        if predictions[i][0] > 0.5:
            predictions[i][0] == 1.0
        else:
            predictions[i][0] == 0.0
        print(predictions[i], "expected ", y[i])
    accuracy = accuracy_score(y, predictions)
    print('Accuracy ((TP + TN) / (P + N)): %f' % accuracy)
    precision = precision_score(y, predictions)
    print('Precision (TP / (TP + FP)): %f' % precision)
    recall = recall_score(y, predictions)
    print('Recall (TP / (TP + FN)): %f' % recall)
    f1 = f1_score(y, predictions)
    print('F1 Score (2TP / (2TP + FP + FN)): %f' % f1)
    predictInput = []
    for kmer in kmerArray:
            encodedKmer = np.append(oneHotEncode(kmer).flatten(), 35)
            predictInput.append(np.array(encodedKmer))
    Xnew = np.array(predictInput)
    predictionsNew = model.predict_classes(Xnew)
    f = open("735.txt", "w")
    for i in range(len(predictionsNew)):
        if predictionsNew[i] == 1:
            f.write(encodeReverse(Xnew[i])+"\n")'''


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trains a model with additional k-mer sets (.txt files) as input.')
    parser.add_argument('-k', metavar='k', type=int, help='K-mer size (k) for the UHS')
    parser.add_argument('--file', metavar= 'file', help='Additional k-mer set file (use preproc.py to merge k-mer sets into one file)')
    parser.add_argument('-e', metavar='e', type=int, help='Number of epochs')
    parser.add_argument('-b', metavar='b', type=int, help='Batch size')

    args = parser.parse_args()
    k = args.k
    UHSfile = args.file
    epochs = args.e
    batchSize = args.b
    additionalPredict(k, UHSfile, epochs, batchSize)
    
    
    
    
    
