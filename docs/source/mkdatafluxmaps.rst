Summary
=======
This app is aimed to prepare the data sample to the Anisotropy analysis. 
Several steps are computed:
   
    Sum in time: 
        the data selection is preferentially done in 1-year-wide time
        bins, so now we want to merge those time bins to get
        the total counts and exposures. Counts and exposure maps
        in each micro energy bin are saved in output/output_counts/.
        Those maps are automatically retrieved if already present in
        the folder and overwrite keyword is False.
    Foreground Subtraction
        If the associated flag is activated the fit and subtraction of the 
        gamma-ray diffuse galactic emission model is performed over the micro
        energy bins. The fit is performed after downgrading the order of the maps 
        and in the count space.
    Flux computation
        The flux maps are derived as the ration of the count maps and the exposure maps
        in each micro energy bin. the flux maps are given in ph/cm2/s/sr. The maps are 
        saved in output/output_flux/.
    Sum in energy
        Flux maps in the macro energy bins are obtained summing the flux maps 
        in the micro bins.
    Save the output file:
       A txt file is automatically generated and saved in output/. It contains 
       all the main parameters related to the run.

.. argparse::
   :filename: ../bin/mkdatafluxmaps.py
   :func: PARSER
   :prog: mkdatafluxmaps.py

   
