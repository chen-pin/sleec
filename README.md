# SLEEC:  Super Lyot ExoEarth Coronagraph

This coding project develops a software suite for designing SLEEC.  SLEEC aims 
to design, implement, and validate the performance of new coronagraphs 
for the (conceptual) HabEx and LUVOIR space observatories.  The Lyot 
coronagraph is SLEEC's foundation, and we intend to advance coronagraphic 
performance by optimizing the complex transmission functions of all the major 
optical planes:  entrance pupil/DM1, DM2, focal, and Lyot pupil planes.  

## Project Status  

### To Do List  
* Write plotting and other utility functions for viewing EMCEE results  

## Getting Started  
Stay tuned.

### Prerequisites  
* Python 3
* Numpy
* Matplotlib.pyplot
* Scipy.sparse
* Scipy.interpolate
* Astropy,io.fits
* HCIPy v0.3.0 w/ bug fix #55 or HCIPy v0.3.1


### Installing
* Install prerequisites (see above)  
* Clone this repository  
* Recommended: 
    * From the code directory <sleec/source>, enter the following commands ("$" represents the command prompt, do not type it):   
        * $ mkdir ../../sleec_out
        * $ mkdir ../../sleec_in

## Running the tests
* From the code directory <sleec/source>, enter the following commands:
    * $ python sleec_opt.py
    * $ python sleec_input.py

### Break down into end to end tests
Stay tuned.

