import numpy as np
import pandas as pd
import itertools
import os
import argparse

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Converts universal hitting sets for given k-mer length intervals into integer encodings (k0, k1).')
	parser.add_argument('-k0', metavar='k0', type=int, help='Minimum k-mer size for the UHS')
	parser.add_argument('-k1', metavar= 'k1', type=int, help='Maximum k-mer size for the UHS')

args = parser.parse_args()
mink=args.k0
maxk=args.k1

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
    if (os.path.isfile("uhs/PASHA" +str(k) + "_" + str(L) + ".txt")):
      f = open("uhs/PASHA" +str(k) + "_" +  str(L) + ".txt", "r")
      g = open("int/PASHA" +str(k) + str(L) +".int", "w")
      for line in f:
        kmer = line.strip("\n")
        totval = 0
        for pos,char in enumerate(kmer):
          val = basedict[char] * pow(4, (k - pos - 1))
          totval += val
        g.write(str(totval) + "\n")
      f.close()
      g.close()



        
