# SLEEC:  Super Lyot ExoEarth Coronagraph

This coding project develops a software suite for designing SLEEC.  SLEEC aims 
to design, implement, and validate the performance of new coronagraphs 
for the (conceptual) HabEx and LUVOIR space observatories.  The Lyot 
coronagraph is SLEEC's foundation, and we intend to advance coronagraphic 
performance by optimizing the complex transmission functions of all the major 
optical planes:  entrance pupil/DM1, DM2, focal, and Lyot pupil planes.  

### Prerequisites  
* Python 3
* Numpy
* Matplotlib.pyplot
* Scipy.sparse
* Scipy.interpolate
* Astropy,io.fits
* HCIPy v0.3.1 with bug fix #70


### Installation
* Install prerequisites (see above)  
* Clone this repository  
* Recommended: 
    * From the code directory <sleec/sleec>, enter the following commands ("$" represents the command prompt, do not type it):  
```  
mkdir ../../sleec_out  
mkdir ../../sleec_in  
```  

#### How to Install HCIPy Bug-Fix Patches  
These instructions assume that one uses the Anaconda Python distribution and 
the Anaconda code directory path is <~/Anaconda3/>.  

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
    * Replace the <hcipy> directory under <~/Anaconda3/> with the downloaded folder    

## Running the tests
* From the code directory <sleec/source>, enter the following command:
    * To test the propagation code:  `python sleec_input.py`.  This does 
not run the optimization routine.  


