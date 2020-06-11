# uhsNN: Generating universal hitting sets using deep neural networks

First, use `preproc.py` to download and preprocess sets:

`python preproc.py -k [k-mer size]`

This will create a `.uhs` file in your directory to train the NN with.

Train the NN and evaluate accuracy with:

`python train.py -k [k-mer size] -f [UHS file] -e [number of epochs] -b [batch size] -o [output model name]`

Example command: `python train.py -k 8 -f DOCKS8_preproc.uhs -e 100 -b 256 -o myModel`

This will create a `.model` file in the output directory.

Then in the main directory, run
`cmake .`<br>
`cd build`<br>
`cmake --build .`

which will compile the prediction executable.

Run the prediction executable with the command

`./predict [model file] [k-mer size] [sequence size] [decycling set output] [additional set output]`

Example command: `./predict myModel.model 8 25 decyc8.txt add825.txt`



