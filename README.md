    
          ---------------------------------------------------
          |  Polaron Coherence in Conjugated Polymer Films  |
          ---------------------------------------------------

What is it?
-----------

This program calculates the Charge Modulation Spectrum in the near IR and the mid-IR region for a 2 Dimensional molecular aggregate by solving the Holstein molecular crystal Hamiltonian using a multiparticle basis set. For elaborate details regarding the theory please go through the following papers:  




How do I use it?
----------------

The program can be compiled using make and the supplied makefile.
You may need to edit the makefile so that it is appropriate for
your system. 

Once compiled, the program can be run from the command line. An
input file containing simulation parameters can be provided as
the first argument. An example input file with parameter explinations
is provided in example.inp. For example, run the program with 

./exciton1D ./example.inp

After running the program, four csv files are produced, one contains
the simulation parameters and has the extension _para.csv while the
other three contain simulation data. The first has the
extension _ab.csv, the second has the extension _disp.csv and the
third has the extension _mom.csv. The first contains the absorption 
spectrum, the second the dispersion curve of the lowest energy exciton,
and the third, the moments of the absorption spectrum. A simple
python script called showspec.py is provided to view the program output. 
To use this script enter the command

python ./showspec.py ./task_title

Where task_title is the name given to the simulation. For example,
in the example.inp file task_title is set to example.


What libraries are required?
----------------------------

To build the exciton1D program:
LAPACK
    
    
    
    Charge Modulation Spectra in 2 dimensional aggregates


    This program calculates the Charge Modulation Spectra for a 2 Dimensional polymer
    aggregate. It is able to run over multiple processors using MPI, Slepc and Petsc, allowing
    for larger system sizes in less time. 
    
    The makefile is attached for compilation.
    

    
    

