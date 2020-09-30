import numpy as np
import pandas as pd
import itertools
import os
mink=10
maxk=10
basedict={"A": 0, "C": 1, "G": 2, "T": 3}
for k in range (mink,maxk+1,1):
  print(k)
  f = open("uhs/decyc" +str(k)+".txt", "r")
  g = open("int/decyc" + str(k) +".int", "w")
  for line in f:
    kmer = line.strip("\n")
    totval = 0
    for pos,char in enumerate(kmer):
      val = basedict[char] * pow(4, (k - pos - 1))
      totval += val
    g.write(str(totval) + "\n")
  f.close()
  g.close()
  for L in range(20,201,10):
    print(L)
    if (os.path.isfile("uhs/DOCKS" +str(k) + "_" + str(L) + ".txt")):
      f = open("uhs/DOCKS" +str(k) + "_" + str(L) + ".txt", "r")
      g = open("int/PDOCKS" +str(k) + str(L) +".int", "w")
      for line in f:
        kmer = line.strip("\n")
        totval = 0
        for pos,char in enumerate(kmer):
          val = basedict[char] * pow(4, (k - pos - 1))
          totval += val
        g.write(str(totval) + "\n")
      f.close()
      g.close()



        
