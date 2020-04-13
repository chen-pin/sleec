# SLEEC:  Super Lyot ExoEarth Coronagraph

This coding project develops a software suite for designing SLEEC.  SLEEC aims 
to design, implement, and validate the performance of new coronagraphs 
for the (conceptual) HabEx and LUVOIR space observatories.  The Lyot 
coronagraph is SLEEC's foundation, and we intend to advance coronagraphic 
performance by optimizing the complex transmission functions of all the major 
optical planes:  entrance pupil/DM1, DM2, focal, and Lyot pupil planes.  

## Project Status  
### To Do List  
1. Integrate update functions from <sleec_input.py> into <sleec_optics.py>
`Coronagraph` class.  
2. Write a method to parse the yml file in <sleec_input.py> and clean up
the rest of <sleec_input.py>.  
3. In <sleec_optics.py>, import/write methods to compute  
    a. single annular contrast  
    b. mean contrast in dark zone  
    c. thruput  
4. In <sleec_input.py>, write cost function  
5. In <sleec_input.py>, write function to do <emcee> optimization.  
6. Add run-logging option    

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

