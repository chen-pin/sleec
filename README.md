# SLEEC:  Super Lyot ExoEarth Coronagraph

This coding project develops a software suite for designing SLEEC.  SLEEC aims 
to design, implement, and validate the performance of new coronagraphs 
for the (conceptual) HabEx and LUVOIR space observatories.  The Lyot 
coronagraph is SLEEC's foundation, and we intend to advance coronagraphic 
performance by optimizing the complex transmission functions of all the major 
optical planes:  entrance pupil/DM1, DM2, focal, and Lyot pupil planes.  

## Prerequisites  
* Python 3
* Numpy
* Matplotlib.pyplot
* Scipy.sparse
* Scipy.interpolate
* Astropy,io.fits
* HCIPy v0.3.1 with bug fix #70


## Installation
* Install prerequisites (see above)  
* Clone this repository  
* Recommended: 
    * From the code directory 'sleec/sleec', enter the following commands ("$" represents the command prompt, do not type it):  
```  
mkdir ../../sleec_out  
mkdir ../../sleec_in  
```  

### How to Install HCIPy Bug-Fix Patches  
These instructions assume that one uses the Anaconda Python distribution and 
the Anaconda code directory path is '~/Anaconda3/'.  

* Via Github cloning  
```  
pip uninstall hcipy  
git clone https://github.com/ehpor/hcipy  
cd hcipy  
pip install -e .  
```  
* By downloading from the Github bug-fix branch
    * Dowload from the zip archive from the branch  
    * Unzip the archive
    * Locate the 'hcipy' directory under '~/Anaconda3/'.  The following are some path-location examples: 
        * /home/pinchen/anaconda3/lib/python3.7/site-packages/hcipy
        * /home/pinchen/anaconda3/envs/dev/lib/python3.8/site-packages/hcipy, if you have an Anaconda environment named "dev."
    * Replace the 'hcipy' directory under '~/Anaconda3/' with the downloaded folder    

## Running tests
* From the code directory <sleec/source>, enter the following command:
    * To test the propagation code:  `python sleec_input.py`.  This does 
not run the optimization routine.  
* Each module has various test at the bottom of the code, some are obsolete.  

## Instructions for adding a coronagraph optimization feature  
* Add function in <sleec_alter.py>.  This is the function that accepts optimization parameters and alters the coronagraph accordingly.  
* In <sleec_pars.yml>  
    * Add `xxx_type` parameter  
    * Add `num_xxx_pars` parameter  
    * Add row(s) in `bound_pars`  
    * Update `nwalkers` value  
    * Update `ndim` value
    * If applicable, add initialization parameter(s)
* In <sleec_run.py>  
    * Add "if" block in <ln_probability()>  
    * Add "if" block in <run_emcee()>  
    * If applicable, add code in <sleec_input.py.import_pars()> to implement new DM feature  
    * If applicable, add code in <sleec_run.ln_prior()> to implement corresponding prior distribution.  
    
Commit <46f32b23> is an example of changes made to implement a new basis-set feature for adjusting DM1.  Use the following command to see the code changes:  git diff 46f32 0bd2e  

# Legal Notices  

Copyright 2025, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.
 
This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons.


