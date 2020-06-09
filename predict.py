import argparse
import math
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Constructs a UHS using predictions from the trained model.')
    parser.add_argument('-k', metavar='k', type=int, help='K-mer size (k) for the UHS')
    parser.add_argument('-L', metavar='L', type=int, help='Sequence size (L) for the UHS')
    parser.add_argument('-m', metavar= 'model', help='Model architecture as JSON file')
    parser.add_argument('-w', metavar= 'weight', help='Model weights as h5 file')


