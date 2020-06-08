import argparse
from requests import get 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Auxiliary code to generate one UHS file for training.')
    parser.add_argument('-k', metavar='k', type=int, help='K-mer size (k) for the UHS')
    parser.add_argument('-m', metavar='m', type=int, help='UHS method - choose from DOCKS (1), DOCKSany (2), and PASHA (3)')

    args = parser.parse_args()
    method = 'DOCKS'
    if args.m == 1:
        method = 'DOCKS'
    elif args.m == 2:
        method = 'DOCKSany'
    elif args.m == 3:
        method = 'PASHA'

    filename = method + str(args.k) + "_preproc.uhs"
    if args.k == 5:
        end = 50
    elif args.k == 6:
        end = 80
    elif args.k == 7:
        end = 120
    elif args.k == 8:
        end = 150
    else:
        end = 210
    with open(filename, "w") as file:
        file.write(">" + 'decycling' + '\n')
        file.close()
    with open(filename, "ab") as file:
            url = 'http://pasha.csail.mit.edu/sets_july24/k' + str(args.k) + '/' + 'decyc' + str(args.k) + '.txt'
            print('Reading from ' + url)
            # get request
            response = get(url)
            # write to file
            file.write(response.content)
            file.close()
    for l in range(20, end, 10):
        with open(filename, "a") as file:
            file.write(">" + str(l) + '\n')
        file.close()
        with open(filename, "ab") as file:
            url = 'http://pasha.csail.mit.edu/sets_july24/k' + str(args.k) + '/' + method + str(args.k) + '_' + str(l) + '.txt'
            print('Reading from ' + url)
            # get request
            response = get(url)
            # write to file
            file.write(response.content)
            file.close()