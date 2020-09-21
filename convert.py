import numpy as np
import pandas as pd
import itertools
import os
mink=5
maxk=10
basedict={"A": 0, "C": 1, "G": 2, "T": 3}
for k in range (mink,maxk+1,1):
  print(k)
  f = open("uhs/k" + str(k) + "/decyc" +str(k)+".txt", "r")
  g = open("decyc" + str(k) +".int", "w")
  for line in f:
    kmer = line.strip("\n")
    totval = 0
    for pos,char in enumerate(kmer):
      val = basedict[char] * pow(4, (k - pos - 1))
      totval += val
    g.write(str(totval) + "\n")
  f.close()
  g.close()
  for L in range(20,210,10):
    print(L)
    if (os.path.isfile("uhs/k" + str(k) + "/DOCKS" +str(k)+ "_" + str(L) + ".txt")):
      f = open("uhs/k" + str(k) + "/DOCKS" +str(k)+ "_" + str(L) + ".txt", "r")
      g = open("DOCKS" +str(k)+ "_" + str(L) +".int", "w")
      for line in f:
        kmer = line.strip("\n")
        totval = 0
        for pos,char in enumerate(kmer):
          val = basedict[char] * pow(4, (k - pos - 1))
          totval += val
        g.write(str(totval) + "\n")
      f.close()
      g.close()



        
