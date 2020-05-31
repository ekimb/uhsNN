# uhsNN: Generating universal hitting sets using deep neural networks

First, generate a decycling set with `pasha`:

`./pasha 7 20 24 decyc7.txt hit720.txt`

Use decycling set to train the NN and evaluate accuracy with:

`python uhsNN.py decyc7.txt 7`

Resulting in output

`Accuracy ((TP + TN) / (P + N)): 0.996399
Precision (TP / (TP + FP)): 0.989713
Recall (TP / (TP + FN)): 0.985068
F1 Score (2TP / (2TP + FP + FN)): 0.987385`

Decycling sets for k = 7 to k = 11 are provided.



