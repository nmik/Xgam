Installation and setup
======================


Required softwares
------------------

* Install the Fermitools (https://github.com/fermi-lat/Fermitools-conda)
* Install PolSpice
* Install healpy python package


Setting up the environment
--------------------------

Set the environment in your .bashrc or .bash_profile:

    function Xgam_env {
    	echo 'Setting Xgam environment...'
    	conda activate fermi
    	export PYTHONPATH=:/path/to/this/package/:${PYTHONPATH}
   	 	export PATH=/path/to/this/package/Xgam/bin:${PATH}
    	export P8_DATA=/path/to/data_files
    	echo 'done.'}

Change the directory where the data files are stored in __init__.py (line 35):

    FT_DATA_FOLDER = '/path/to/data_files'

In the data directory should there be the following folders:
   
   * photon/      -> where FT1 files are stored
   * spacecraft/  -> where FT2 files are stored
   * output/      -> where ST outputs will be stored


Docker alternative
------------------
See the dedicated wiki page.

How to run the analysis
-----------------------
This will show all the possible settings of a given function in bin/ :

     python bin/mkAnyFunction.py -h 

Some more details can be found in the Wiki of this project:

https://github.com/nmik/Xgam/wiki

(A more detailed documentation will be soon available).
