# uhsNN: Generating universal hitting sets using deep neural networks

First, generate a decycling set with `pasha`:

`./pasha 7 20 24 decyc7.txt hit720.txt`

Use decycling set to train the NN and evaluate accuracy with:

`python uhsNN.py decyc7.txt 7`

Resulting in output

`Accuracy: 99.82`




