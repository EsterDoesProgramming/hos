# HOS
High order spectral (HOS) code for 3D water wave simulations.

Contents:

* Csource/ Source files and XCode project
* matlab/ Matlab utilities for pre- and post-processing
* python/ Python utilities for pre- and post-processing

This software uses two external libraries:

* FFTW (www.fftw.org)
* HDF5 (http://www.hdfgroup.org/)

## Getting started

Get a copy of the code using `git clone`.
Several test cases have been implemented to demonstrate the model's functionalities.
Before you start working with the code, it is advisable to run those and look at the results.
There are currently two ways of doing this: A. Using IPython notebooks; and B. Using the terminal alone.

But first, please configure the Makefile at `~/hos/Csource/2dpar` for the parallel version to fit your local system.  In particular, check that the library paths are set correctly. Hint: The terminal command `locate` is your friend in case you are unsure about where your libraries are. 

### A. The IPython Notebook Way
Go to `~/hos/python` and start your ipython notebook server. If you are new to IPython notebooks `https://ipython.org/notebook.html` offers comprehensive information on how to do that.
The ipython notebook `HOS_Testsuite.ipynb` has a documentation of all test cases and can be used to compile and run the code.
Just follow the instructions given in the notebook.

### B. Using Terminal
If you are not a python fan there are two things to do: A. It is strongly recommended that you reconsider your opinion :); and B. You can find a `Makefile.Testcase` in `~/hos/Csource/2dpar/` that has the set up for all the testcases. Just use `make <NameOfTestCase>` to compile and run. 

### C. Create Your Own Initial Conditions
To create your own initial conditions, go to `~/hos/python` and start your ipython notebook server. The notebook `HOS_InitialCondition.ipynb` can be used to define the initial conditions.
Several commonly used, analytical spectra have already been implemented and documented. However, feel free to add your own!
The HOS model currently supports two forms of input parameters: ASCII and HDF5. Both formats will be written by the notebook.


The code is documented with doxygen.
To obtain the full version of the documentation for the parallel version of the code (in html and LaTeX format), simply open a terminal, go to ~/hos/Csource/2dpar/ and run

> doxywizard HOSM.doxy

This creates a folder _documentation_, that contains all the information you need.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

