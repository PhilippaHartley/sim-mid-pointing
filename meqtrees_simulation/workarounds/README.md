For a workaround to generate a set of pointing errors random in time and across antennae:

* replace  ``` ErrorGens.py``` in meqtrees directory ```Cattery/Siamese/OMS/``` (or just add the function ```Random_Dynamic``` from this version of the file).

For a workaround to write out a numpy array -- ```meqtrees_dlm.npy``` -- of pointing errors with the shape (2, mxn), where each column contains 2 pointing errors, dl and dm, and each row contains values ordered first by by timestamp m and then by antenna m:

* replace ```__init__.py``` in timba directory ```Grid/```. For this to work, the simulation must be run from the GUI,  a bookmark for ``E ponting error: by station`` must be opened after compilation and before simulation run, and the tile size (in the run simulation window) set to 65 (in order to match the current number of timestamps at the time of writing). There is probably a much more elegant way of doing this using PyNodes.
