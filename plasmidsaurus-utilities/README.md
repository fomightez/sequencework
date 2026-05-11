Plasmidsaurus-utilities
-----------------------

My collection of items to help with sequencing using the Plasmidsaurus service.

Dealing with things when you get your 'Huge Plasmid'/BAC sequences back from Plasmidsaurus
--------------------------

Plasmidsaurus' related info is at ['Whole Genome Technical Documentation'](https://plasmidsaurus.com/technical-documentation/genome)

related discussion I found after I started my own plans, works for AmpR plasmids only it seems:
['Plasmidsaurus fasta standardizer' Posted on April 6, 2023 by kmatreyek](https://www.matreyeklab.com/plasmidsaurus-fasta-standardizer/1803/)

See my notebook `Rotate SLASH Permute the sequences to match your favorite BAC.ipynb`.  
The idea is that you work through that as it is and then next modify it to use yourself with your own sequences in the same or a new session.  
Presently, you can run that notebook in MyBinder-served sessions that include biopython installed.  
I have, though, set up things so that you can launch using the following launch badge and not only will the environment be set up appropriately, all the necessary ancillary files will also be present:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fomightez/sequencework/HEAD?urlpath=%2Flab%2Ftree%2Fplasmidsaurus-utilities%2FRotate%20SLASH%20Permute%20the%20sequences%20to%20match%20your%20favorite%20Huge%20Plasmid%20or%20BAC.ipynb)

Please use that badge or [click here to directly launch the recommended session](https://mybinder.org/v2/gh/fomightez/sequencework/HEAD?urlpath=%2Flab%2Ftree%2Fplasmidsaurus-utilities%2FRotate%20SLASH%20Permute%20the%20sequences%20to%20match%20your%20favorite%20Huge%20Plasmid%20or%20BAC.ipynb) to get things going in the easiest way possible.

If you are an advanced user and just want a more basic session with biopython and are comfortbale having to get & place the ancillary in the running session yourself, then here are some additional options where Biopython is installed and working:  
- [my cl_demo-binder) here](https://github.com/fomightez/cl_demo-binder)
- (I'm not sure launches from the seemingly more appropriate [repo 'cl_sq_demo-binder'](https://github.com/fomightez/cl_sq_demo-binder) are presently 
failing.)

Jupyter interfaces with Biopython already installed:   
JupyterLab interface: [![Binder](https://mybinder.org/v2/gh/fomightez/sequencework/HEAD?urlpath=%2Flab%2Ftree%2Fplasmidsaurus-utilities%2FRotate%20SLASH%20Permute%20the%20sequences%20to%20match%20your%20favorite%20Huge%20Plasmid%20or%20BAC.ipynb)  
Jupyter Notebook 7+:  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fomightez/sequencework/HEAD?urlpath=%2Ftree%2Fplasmidsaurus-utilities%2FRotate%20SLASH%20Permute%20the%20sequences%20to%20match%20your%20favorite%20Huge%20Plasmid%20or%20BAC.ipynb)

**(Still to explore: Is JupyterLite also an option for running the demo notebook?)**
