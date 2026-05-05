Plasmidsaurus-utilities
-----------------------

My collection of items to help with sequencing using the Plasmidsaurus service.

Dealing with things when you get your BAC sequences back from Plasmidsaurus
--------------------------

Plasmidsaurus' related info is at ['Whole Genome Technical Documentation'](https://plasmidsaurus.com/technical-documentation/genome)

related discussion I found after I started my own plans, works for AmpR plasmids only it seems:
['Plasmidsaurus fasta standardizer' Posted on April 6, 2023 by kmatreyek](https://www.matreyeklab.com/plasmidsaurus-fasta-standardizer/1803/)

See my notebook `Rotate SLASH Permute the sequences to match your favorite BAC.ipynb`.  
The idea is that you work through that as it is and then next modify it to use yourself with your own sequences.
Presently can run it at sessions launched from [here](https://github.com/fomightez/cl_demo-binder) that include biopython installed. [Click here to launch a session from there in the JupyterLab interface]([cl_demo-binder/master](https://mybinder.org/v2/gh/fomightez/cl_demo-binder/master?urlpath=%2Flab%2Ftree%2Findex.ipynb), which is recommended for this notebook since you'll want to drag-and-drop in your sequence file. (I'm not sure launches from the seemingly more appropriate [repo 'cl_sq_demo-binder'](https://github.com/fomightez/cl_sq_demo-binder) are presently failing.)

JupyterLab interface with Biopython already installed:   
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fomightez/cl_demo-binder/master?urlpath=%2Flab%2Ftree%2Findex.ipynb)  

(Still to explore: Is JupyterLite also an option for running the demo?)
