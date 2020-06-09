import argparse
import math
import graph
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Constructs a UHS using predictions from the trained model.')
    parser.add_argument('-k', metavar='k', type=int, help='K-mer size (k) for the UHS')
    parser.add_argument('-L', metavar='L', type=int, help='Sequence size (L) for the UHS')
    parser.add_argument('-m', metavar= 'model', help='Model architecture as JSON file')
    parser.add_argument('-w', metavar= 'weight', help='Model weights as h5 file')
    parser.add_argument('-d', metavar= 'decyc', help='Decycling set as txt file')

    args = parser.parse_args()
    graph = graph.Graph(args.k)
    print("Graph constructed.")
    print('Removing decycling set...')
    graph.removeDecycling(args.d)
    print('Starting additional set...')
    graph.HittingPipeline(args.L)
    print('Done!')
