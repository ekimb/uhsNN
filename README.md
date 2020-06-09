# uhsNN: Generating universal hitting sets using deep neural networks

First, use `preproc.py` to download and preprocess sets:

`python preproc.py -k [k-mer size]`

This will create a `.uhs` file in your directory to train the NN with.

Train the NN and evaluate accuracy with:

`python uhsNN.py -k [k-mer size] -f [UHS file] -e [number of epochs] -b [batch size]`



