    
          ---------------------------------------------------
          |  Polaron Coherence in Conjugated Polymer Films  |
          ---------------------------------------------------

What is it?
-----------

This program calculates the Charge Modulation Spectrum in the near IR and the mid-IR region for a 2 Dimensional molecular aggregate by solving the disordered Holstein molecular crystal Hamiltonian employing a multiparticle basis set. 

Modelling Disorder in Polymer Aggregates
----------------------------------------

Energetic disorder in polymer films arises from variations in nearest neighbor packing distances, change in intramolecular torsional angles and the presence of spatially fluctuating electric fields. We introduce various kinds of diagonal and offdiagonal disorder. Please look into the reference mentioned below for the disorder models that have been used in the code. 


Polaron Absorption and Coherence 
----------------------------------------

We have developed a relation between polaron coherence and polaron absorption spectra. To put it in simple terms, given an absorption spectra we can extract the polaron delocalization lengths both along the polymer backbone and in between chains. The subroutine to calculate the polaron absorption and coherence has been included in the repository. Please have a look

MPI  
----

The program uses MPI techniques to parallely run the disorder configurations.
    

    
    

