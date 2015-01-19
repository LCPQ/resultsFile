# resultsFile
Python interface to read output files of quantum chemistry programs

To add a module to read a new kind of output file, just add a file 
in the `Modules` directory.



# Using the library

Example (`resultsFile` is supposed to be in your `sys.path`):

``` Python
import resultsFile 

file = resultsFile.getFile("g09_output.log")
print 'recognized as', str(file).split('.')[-1].split()[0]
print file.mo_sets

```

## Constraints

### Gaussian09

* `GFPRINT`  : Needed to read the AO basis set
* `pop=Full` : Needed to read all the MOs
* `#p CAS(SlaterDet)` : CAS-SCI CI coefficients

When doing a CAS with Gaussian, first do the Hartree-Fock calculation saving the checkpoint
file and then do the CAS in a second calculation.


### Molpro

* `print, basis;` :  Needed to read the AO basis set
* `gprint,orbital;` : Needed to read the MOs
* `gprint,civector; gthresh,printci=0.;` : Needed to read the CI coefficients
* `orbprint` : Ensures all the MOs are printed

An RHF calculation is mandatory before any MCSCF calculation, since some
information is printed only the RHF section. Be sure to print *all* molecular
orbitals using the `orbprint` keyword, and to use the same spin multiplicity
and charge between the RHF and the CAS.


### GAMESS-US

For MCSCF calculations, first compute the MCSCF single-point wave function with
the GUGA algorithm. Then, put the the MCSCF orbitals (of the `.dat` file) in
the GAMESS input file, and run a single-point GUGA CI calculation with the
following keywords:

* `PRTTOL=0.0001` in the `$GUGDIA` group to use a threshold of 1.E-4 on the CI coefficients
* `NPRT=2` in the `$CIDRT` group to print the CSF expansions in terms of Slater determinants
* `PRTMO=.T.` in the `$GUESS` group to print the molecular orbitals

# Debugging

Any module can be run as an stand-alone executable. For example:

```
$ resultsFile/Modules/gamessFile.py

    resultsFile version 1.0, Copyright (C) 2007 Anthony SCEMAMA
    resultsFile comes with ABSOLUTELY NO WARRANTY; for details see the
    gpl-license file.
    This is free software, and you are welcome to redistribute it
    under certain conditions; for details see the gpl-license file.

Usage:
------

resultsFile/Modules/gamessFile.py [options] file

Options:
--------

  --date                      :  When the calculation was performed.
  --version                   :  Version of the code generating the file.
  --machine                   :  Machine where the calculation was run.
  --memory                    :  Requested memory for the calculation.
  --disk                      :  Requested disk space for the calculation.
  --cpu_time                  :  CPU time.
  --author                    :  Who ran the calculation.
  --title                     :  Title of the run.
  --units                     :  Units for the geometry (au or angstroms).
  --methods                   :  List of calculation methods.
  --options                   :  Options given in the input file.
  --spin_restrict             :  Open-shell or closed-shell calculations.
  --conv_threshs              :  List of convergence thresholds.
  --energies                  :  List of energies.
  --one_e_energies            :  List of one electron energies.
  --two_e_energies            :  List of two electron energies.
  --ee_pot_energies           :  List of electron-electron potential energies.
  --Ne_pot_energies           :  List of nucleus-electron potential energies.
  --pot_energies              :  List of potential energies.
  --kin_energies              :  List of kinetic energies.
  --virials                   :  Virial ratios.
  --point_group               :  Symmetry used.
  --num_elec                  :  Number of electrons.
  --charge                    :  Charge of the system.
  --multiplicity              :  Spin multiplicity of the system.
  --nuclear_energy            :  Repulsion of the nuclei.
  --dipole                    :  Dipole moment
  --geometry                  :  Atom types and coordinates.
  --basis                     :  Basis set definition
  --mo_sets                   :  List of molecular orbitals
  --mo_types                  :  Types of molecular orbitals (canonical, natural,...)
  --mulliken_mo               :  Mulliken atomic population in each MO.
  --mulliken_ao               :  Mulliken atomic population in each AO.
  --mulliken_atom             :  Mulliken atomic population.
  --lowdin_ao                 :  Lowdin atomic population in each AO.
  --mulliken_atom             :  Mulliken atomic population.
  --lowdin_atom               :  Lowdin atomic population.
  --two_e_int_ao              :  Two electron integrals in AO basis
  --determinants              :  List of Determinants
  --num_alpha                 :  Number of Alpha electrons.
  --num_beta                  :  Number of Beta electrons.
  --closed_mos                :  Closed shell molecular orbitals
  --active_mos                :  Active molecular orbitals
  --virtual_mos               :  Virtual molecular orbitals
  --determinants_mo_type      :  MO type of the determinants
  --det_coefficients          :  Coefficients of the determinants
  --csf_mo_type               :  MO type of the determinants
  --csf_coefficients          :  Coefficients of the CSFs
  --symmetries                :  Irreducible representations
  --occ_num                   :  Occupation numbers
  --csf                       :  List of Configuration State Functions
  --num_states                :  Number of electronic states
  --two_e_int_ao_filename     :  
  --one_e_int_ao_filename     :  
  --atom_to_ao_range          :  
  --gradient_energy           :  Gradient of the Energy wrt nucl coord.
  --text                      :  
  --uncontracted_basis        :  
  --uncontracted_mo_sets      :  

```

