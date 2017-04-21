#!/usr/bin/python
#   resultsFile is a library which allows to read output files of quantum 
#   chemistry codes and write input files.
#   Copyright (C) 2007 Anthony SCEMAMA 
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program; if not, write to the Free Software Foundation, Inc.,
#   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#   Anthony Scemama
#   LCPQ - IRSAMC        
#   Universite Paul Sabatier
#   118, route de Narbonne      
#   31062 Toulouse Cedex 4      
#   scemama@irsamc.ups-tlse.fr 


from lib import *
import lib.basis as Basis
import sys
import copy

fileTypes = []

local_vars = [ \
   # File properties
    ( 'filename'     , "Name of the results file."),
    ( 'text'         , "Text of the results file."),
    ( 'pos'          , "Position in the results file."),
    ( 'date'         , "When the calculation was performed."),
    ( 'version'      , "Version of the code generating the file."),
    ( 'author'       , "Who ran the calculation."),
    ( 'machine'      , "Machine where the calculation was run."),
    ( 'memory'       , "Requested memory for the calculation."),
    ( 'disk'         , "Requested disk space for the calculation."),
    ( 'num_proc'     , "Number of processors used."),
    ( 'cpu_time'     , "CPU time."),
   # General properties
    ( 'title'         , "Title of the run."),
    ( 'options'       , "Options given in the input file."),
    ( 'units'         , "Units for the geometry (au or angstroms)."),
    ( 'methods'       , "List of calculation methods."),
    ( 'spin_restrict' , "Open-shell or closed-shell calculations."),
    ( 'conv_threshs'  , "List of convergence thresholds."),
    ( 'energies'      , "List of energies."),
    ( 'one_e_energies', "List of one electron energies."),
    ( 'two_e_energies', "List of two electron energies."),
    ( 'ee_pot_energies',"List of electron-electron potential energies."),
    ( 'Ne_pot_energies',"List of nucleus-electron potential energies."),
    ( 'pot_energies'  , "List of potential energies."),
    ( 'kin_energies'  , "List of kinetic energies."),
    ( 'virials'       , "Virial ratios."),
    ( 'mulliken_mo'   , "Mulliken atomic population in each MO."),
    ( 'mulliken_ao'   , "Mulliken atomic population in each AO."),
    ( 'lowdin_ao'     , "Lowdin atomic population in each AO."),
    ( 'mulliken_atom' , "Mulliken atomic population."),
    ( 'lowdin_atom'   , "Lowdin atomic population."),
    ( 'dipole'        , "Dipole moment"),
    ( 'quadrupole'    , "Quadrupole moment"),
    ( 'num_states'    , "Number of electronic states"),
   # Geometry properties
    ( 'point_group'    , "Symmetry used."),
    ( 'geometry'       , "Atom types and coordinates."),
    ( 'symmetries'     , "Irreducible representations"),
    ( 'num_elec'       , "Number of electrons."),
    ( 'num_alpha'      , "Number of Alpha electrons."),
    ( 'num_beta'       , "Number of Beta electrons."),
    ( 'charge'         , "Charge of the system."),
    ( 'multiplicity'   , "Spin multiplicity of the system."),
    ( 'nuclear_energy ', "Repulsion of the nuclei."),
    ( 'gradient_energy', "Gradient of the Energy wrt nucl coord."),
   # Basis set
    ( 'basis'             , "Basis set definition"),
    ( 'uncontracted_basis', "Uncontracted Basis set"),
   # Pseudopotentials
    ( 'pseudo'        , "Pseudopotential data"),
   # Orbitals
    ( 'mo_sets'       , "List of molecular orbitals"),
    ( 'mo_types'      , "Types of molecular orbitals (canonical, natural,...)"),
    ( 'occ_num'       , "Occupation numbers"),
    ( 'uncontracted_mo_sets', "List of molecular orbitals in the uncontracted basis set."),
   # Determinants
    ( 'closed_mos'  , "Closed shell molecular orbitals"),
    ( 'active_mos'  , "Active molecular orbitals"),
    ( 'virtual_mos'  , "Virtual molecular orbitals"),
    ( 'csf'           , "List of Configuration State Functions"),
    ( 'determinants'  , "List of Determinants"),
    ( 'csf_mo_type'  , "MO type of the determinants"),
    ( 'determinants_mo_type'  , "MO type of the determinants"),
    ( 'csf_coefficients', "Coefficients of the CSFs"),
    ( 'det_coefficients', "Coefficients of the determinants"),
   # Integrals
    ( 'one_e_int_ao'     , "One electron integrals in AO basis"),
    ( 'two_e_int_ao'     , "Two electron integrals in AO basis"),
    ( 'one_e_int_mo'     , "One electron integrals in MO basis"),
    ( 'two_e_int_mo'     , "Two electron integrals in MO basis"),
   ]

resultsFile_defined_vars = ["text","uncontracted_basis", "uncontracted_mo_sets"]

class resultsFile(object):
   """ Class containing the definition of files.
   """

   local_vars = list(local_vars)
   defined_vars = list(resultsFile_defined_vars)

   def __init__(self,name):
       """All local variables are set to None, except filename.
       """
       for var, doc in local_vars:
          setattr(self,'_'+var,None)
       self._filename = name


   def get_text(self):
       """Reads the whole file.
       """
       if self._text is None:
          try:
             file = open(self.filename,"r")
             self._text = file.readlines()
          except IOError:
             print "Unable to open "+self.filename
             sys.exit(1)
          file.close()
       return self._text

   def get_uncontracted_basis(self):
     if self._uncontracted_basis is None:
       try:
         from resultsFile_cython import get_uncontracted_basis as f
         has_cython = True
       except ImportError:
         has_cython = False
       basis = self.basis
       if has_cython:
         self._uncontracted_basis = f(basis)
       else:
         uncontr = []
         for contr in basis:
            for b in contr.prim:
               uncontr.append(b)
         self._uncontracted_basis = uncontr
     return self._uncontracted_basis

   def get_uncontracted_mo_sets(self):
     if self._uncontracted_mo_sets is None:
       try:
         from resultsFile_cython import get_uncontracted_mo_sets as f
         has_cython = True
       except ImportError:
         has_cython = False
       if has_cython:
         self._uncontracted_mo_sets = f(self.basis,
                         self.uncontracted_basis,
                         self.mo_sets,
                         self.mo_types)
       else:
         uncontr = {}
         basis = self.basis
         for motype in self.mo_types:
           uncontr[motype] = []
           for mo in self.mo_sets[motype]:
             lenmovector = len(mo.vector)
             monew = orbital()
             monew.basis = self.uncontracted_basis
             monew.eigenvalue = mo.eigenvalue
             monew.set = motype
             v = []
             for i, contr in enumerate(basis):
               if i<lenmovector:
                 ci = mo.vector[i]
                 if ci == 0.:
                   for p in range(len(contr.prim)):
                      v.append(0.)
                 else:
                   for p, c in zip(contr.prim,contr.coef):
                      v.append(c*ci/p.norm)
             monew.vector = v
             uncontr[motype].append(monew)
         self._uncontracted_mo_sets = uncontr
       #endif
     return self._uncontracted_mo_sets

   def clean_contractions(self):
     basis = self.basis
     newbasis = []
     idx = range(len(basis))
     for k,b1 in enumerate(basis):
       addBasis=True
       for l, b2 in enumerate(basis[:k]):
         if b2 == b1:
           idx[k] = l
           addBasis=False
           break
       if addBasis:
         newbasis.append(b1)
     self._basis = newbasis 

     mo_sets = self.mo_sets
     for motype in self.mo_types:
       for mo in mo_sets[motype]:
         lenmovector = len(mo.vector)
         newvec = [None for i in idx]
         for i in idx:
           newvec[i] = 0.
         for k,l in enumerate(idx):
           if k < lenmovector:
             newvec[l] += mo.vector[k]
         mo.vector = []
         for c in newvec:
           if c is not None:
             mo.vector.append(c)

   def clean_uncontractions(self):
     basis = self.uncontracted_basis
     newbasis = []
     idx = range(len(basis))
     for k,b1 in enumerate(basis):
       addBasis=True
       for l, b2 in enumerate(basis[:k]):
         if b2 == b1:
           idx[k] = l
           addBasis=False
           break
       if addBasis:
         newbasis.append(b1)
     self._uncontracted_basis = newbasis 

     mo_sets = self.uncontracted_mo_sets
     for motype in self.mo_types:
       for mo in mo_sets[motype]:
         lenmovector = len(mo.vector)
         newvec = [None for i in idx]
         for i in idx:
           newvec[i] = 0.
         for k,l in enumerate(idx):
           if k < lenmovector:
             newvec[l] += mo.vector[k]
         mo.vector = []
         for c in newvec:
           if c is not None:
             mo.vector.append(c)
       
   def convert_to_cartesian(self):
     basis = self.basis
     newbasis = []
     idx = range(len(basis))
     map = []
     weight = []
     for i,b in enumerate(basis):
       l, m = Basis.get_lm(b.sym)
       if l is None:
         newbasis.append(b)
         map.append(i)
         weight.append(1.)
       else:
         powers, coefs = xyz_from_lm(l,m)
         for j,prim in enumerate(b.prim):
           b.coef[j] /= prim.norm
         for c, p in zip(coefs, powers):
           contr = copy.deepcopy(b)
           sym = ''
           for l,letter in enumerate('xyz'):
             sym += p[l]*letter
           contr.sym = sym
           for j,prim in enumerate(contr.prim):
             prim.sym = sym
             contr.coef[j] *= prim.norm
           newbasis.append(contr)
           map.append(i)
           weight.append(c)

     mo_sets = self.mo_sets
     for motype in self.mo_types:
       for mo in mo_sets[motype]:
         newvec = []
         vec = mo.vector
         for i,w in zip(map,weight):
           newvec.append(vec[i]*w)
         mo.vector = newvec

     same_as = {}
     for i,b1 in enumerate(newbasis):
       for j,b2 in enumerate(newbasis[:i]):
         if b1 == b2:
          same_as[i] = j
          weight[j] += weight[i]
          break
     to_remove = same_as.keys()
     to_remove.sort()
     to_remove.reverse()
     for i in to_remove:
       newbasis.pop(i)
       weight.pop(i)
       map.pop(i)


     for motype in self.mo_types:
       for mo in mo_sets[motype]:
         for i in to_remove:
           index = same_as[i]
           value = mo.vector.pop(i)
           mo.vector[index] += value

     self._basis = newbasis
     self._mo_sets = mo_sets

   def find_string(self,chars):
       """Finds the 1st occurence of chars.
       """
       self._pos = 0
       self.find_next_string(chars)

   def find_last_string(self,chars):
       """Finds the 1st occurence of chars.
       """
       self._pos = len(self.text)-1
       self.find_prev_string(chars)

   def find_next_string(self,chars):
       """Finds the next occurence of chars.
       """
       pos  = self._pos
       text = self.text
       found = False
       while not found and pos < len(text):
          if chars in text[pos]:
             found = True
          else:
             pos += 1
       if not found:
          raise IndexError
       self._pos = pos

   def find_prev_string(self,chars):
       """Finds the next occurence of chars.
       """
       pos  = self._pos
       text = self.text
       found = False
       while not found and pos < len(text):
          if chars in text[pos]:
             found = True
          else:
             pos -= 1
       if not found:
          raise IndexError
       self._pos = pos


   for i, j in local_vars:
      if i not in defined_vars:
         exec build_get_funcs(i) in locals()
      exec build_property(i,j) in locals()
   del i,j


def main(fileType):
   import getopt
   print """
    resultsFile version 1.0, Copyright (C) 2007 Anthony SCEMAMA
    resultsFile comes with ABSOLUTELY NO WARRANTY; for details see the
    gpl-license file.
    This is free software, and you are welcome to redistribute it
    under certain conditions; for details see the gpl-license file."""

   full_list = fileType.defined_vars + resultsFile.defined_vars
   try:
     opts, args = getopt.gnu_getopt(sys.argv[1:],'',full_list)
   except getopt.GetoptError:
     args = []
   if len(args) == 0:
     usage(fileType)
     sys.exit(2)
   f = fileType(args[0])
   for o,a in opts:
     print o[2:]
     print ''.join(['-' for k in o[2:]])
     PrettyPrint = prettyPrint
     exec 'PrettyPrint(f.'+o[2:]+')' in locals()
     print ""
   sys.exit(0)

def usage(fileType):
   print ""
   print "Usage:"
   print "------"
   print ""
   print sys.argv[0], '[options] file'
   print ""
   print "Options:"
   print "--------"
   print ""
   for o in fileType.defined_vars + resultsFile.defined_vars:
     line = ("  --"+o).ljust(30)+':  '
     for l in fileType.local_vars:
       if l[0] == o:
         line += l[1]
         break
     print line
   print ""

   
if __name__ == '__main__':
   main(resultsFile)

