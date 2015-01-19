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



import include
eval(include.code)

import struct
import re

import os
QMCCHEM_PATH = os.getenv("QMCCHEM_PATH",default="")
if QMCCHEM_PATH == "":
 print "QmcChem new files are not handled."
 class qmcchem_newFile(resultsFile):
   pass
else:
 sys.path = [ QMCCHEM_PATH+"/scripts" ]+sys.path
 from ezfio import ezfio

 qmcchem_newFile_defined_vars = [ "date", "version", \
                "title", "units", "methods", \
                "point_group", "num_elec", \
                "charge", "multiplicity","geometry",\
                "basis","mo_sets","mo_types",\
                "determinants", "num_alpha", "num_beta",\
                "closed_mos", "active_mos", "virtual_mos", \
                "determinants_mo_type", "det_coefficients", \
                "csf_mo_type", "csf_coefficients", "occ_num", \
                "csf" ]

 class qmcchem_newFile(resultsFile):
   """ Class defining the qmcchem_new file.
   """

   local_vars = list(local_vars)
   defined_vars = list(qmcchem_newFile_defined_vars)

   def __init__(self,name):
     resultsFile.__init__(self,name)
     ezfio.set_filename(self.filename)

   def get_version(self):
     if self._version is None:
       self._version = ezfio.get_version()
     return self._version

   def get_date(self):
     if self._date is None:
       self._date = ezfio.get_ezfio_creation()
     return self._date

   def get_num_elec(self):
     if self._num_elec is None:
       self._num_elec = self.num_alpha + self.num_beta
     return self._num_elec

   def get_multiplicity(self):
     if self._multiplicity is None:
       self._multiplicity = self.num_alpha - self.num_beta + 1
     return self._multiplicity

   def get_charge(self):
     if self._charge is None:
       self._charge = sum(ezfio.get_nuclei_nucl_charge())-float(self.num_elec)
     return self._charge

   def get_title(self):
     if self._title is None:
       self._title = self.filename
     return self._title

   def get_units(self):
     if self._units is None:
       self._units = 'BOHR'
       #for a gamess use the units is give by an option on gamess_write_contrl
     return self._units

   def get_methods(self):
     if self._methods is None:
        self._methods = ['QMC']
     return self._methods

   def get_point_group(self):
     if self._point_group is None:
        self._point_group = "C1"
     return self._point_group 

   def get_geometry(self):
     if self._geometry is None:
       self.get_geometryX()
       self.get_basisX()
     return self._geometry

   def get_basis(self):
     if self._basis is None:
       self.get_geometry()
     return self._basis

   def get_geometryX(self):
     self._geometry = []
     charge = ezfio.get_nuclei_nucl_charge()
     coord  = ezfio.get_nuclei_nucl_coord()
     num    = ezfio.get_nuclei_nucl_num()
     for i in range(num):
        temp = atom()
        temp.charge = charge[i]
        temp.coord  = (coord[0][i], coord[1][i], coord[2][i])
        temp.name   = 'X'
        temp.basis  = []
        self._geometry.append(temp)

   def get_basisX(self):
        coef = ezfio.get_ao_basis_ao_coef()
        expo = ezfio.get_ao_basis_ao_expo()
        nucl = ezfio.get_ao_basis_ao_nucl()
        num  = ezfio.get_ao_basis_ao_num()
        power= ezfio.get_ao_basis_ao_power()
        prim_num= ezfio.get_ao_basis_ao_prim_num()

        self._basis = []
        for i in range(num):
            contr = contraction()
            for j in range(prim_num[i]):
              gauss = gaussian()
              atom = self._geometry[nucl[i]-1]
              gauss.center = atom.coord
              gauss.expo = expo[j][i]
              name = normalize_basis_name('x'*power[0][i]+'y'*power[1][i]+'z'*power[2][i])
              if name == '': name = 's'
              gauss.sym  = name
              contr.append(coef[j][i],gauss)
            self._geometry[nucl[i]-1].basis.append(contr)
            self._basis.append(contr)


   def get_mo_types(self):
      if self._mo_types is None:
        self._mo_types = ['QMC']
      return self._mo_types 

   def get_mo_sets(self):
      if self._mo_sets is None:
         self._mo_sets = {}
         self._mo_sets['QMC'] = []
         coef = ezfio.get_mo_basis_mo_coef()
         energy = ezfio.get_mo_basis_mo_energy()
         num = ezfio.get_mo_basis_mo_tot_num()
         for i in range(num):
           v = orbital()
           v.basis = self.basis
           v.set = 'QMC'
           v.eigenvalue = energy[i]
           v.vector = coef[i]
           self.mo_sets['QMC'].append(v)
      return self._mo_sets
         
   def get_num_alpha(self):
      if self._num_alpha is None:
        self._num_alpha = ezfio.get_electrons_elec_alpha_num()
      return self._num_alpha 

   def get_num_beta(self):
      if self._num_beta is None:
        self._num_beta = ezfio.get_electrons_elec_beta_num()
      return self._num_beta 

   def get_determinants_mo_type(self):
      if self._determinants_mo_type is None:
         self._determinants_mo_type = 'QMC'
      return self._determinants_mo_type

   def get_csf_mo_type(self):
      if self._csf_mo_type is None:
         self._csf_mo_type = 'QMC'
      return self._csf_mo_type

   def get_determinants(self):
      if self._determinants is None:
        determinants = []
        if self.csf is not None:
           for csf in self.csf:
              for new_det in csf.determinants:
                 determinants.append(new_det)
        else:
           pass
        if determinants != []:
          self._determinants_mo_type = self.mo_types[-1]
          self._determinants = determinants
      return self._determinants

   def get_csf(self):
      method = self.methods[0]
      if self._csf is None:
        csf = []
        ncore = ezfio.get_mo_basis_mo_closed_num()
        nact  = ezfio.get_mo_basis_mo_active_num()
        core_a = []
        core_b = []
        for i in range(ncore):
           core_a.append(i)
           core_b.append(i)
        num = ezfio.get_determinants_det_num()
        occ = ezfio.get_determinants_det_occ()
        if occ == []:
          occ = [[[0]],[[0]]]
        for i in range(num):
          this_csf = CSF()
          tempcsf_a = core_a + map(lambda x: x-1, occ[0][i])
          tempcsf_b = core_b + map(lambda x: x-1, occ[1][i])
          this_csf.append(1.,tempcsf_a,tempcsf_b)
          csf.append(this_csf)
        if csf != []:
           self._csf = csf
      return self._csf


   def get_closed_mos(self):
      if self._closed_mos is None:
         cls = ezfio.get_mo_basis_mo_classif()
         self._closed_mos = []
         self._virtual_mos = []
         self._active_mos = []
         for i in range(len(cls)):
           if cls[i] == 'c':
             self._closed_mos.append(i)
           elif cls[i] == 'a':
             self._active_mos.append(i)
           elif cls[i] == 'v':
             self._virtual_mos.append(i)
      return self._closed_mos

   def get_virtual_mos(self):
      if self._virtual_mos is None:
        self.get_closed_mos()
      return self._virtual_mos

   def get_active_mos(self):
      if self._active_mos is None:
        self.get_closed_mos()
      return self._active_mos

   def get_det_coefficients(self):
      if self._det_coefficients is None:
        self._det_coefficients = []
        csf = self.csf
        for state_coef in self.csf_coefficients:
           vector = []
           for i,c in enumerate(state_coef):
              for d in csf[i].coefficients:
                 vector.append(c*d)
           self._det_coefficients.append(vector)
      return self._det_coefficients

   def get_csf_coefficients(self):
      if self._csf_coefficients is None:
           self._csf_coefficients = [ [] ]
           self._csf_coefficients[0] = ezfio.get_determinants_det_coef()
      return self._csf_coefficients

   def get_occ_num(self):
     if self._occ_num is None:
      self._occ_num = {}
      motype = 'QMC'
      self._occ_num[motype] = ezfio.get_mo_basis_mo_occ()
     return self._occ_num 

      
fileTypes.append(qmcchem_newFile)
  
if __name__ == '__main__':
   main(qmcchem_newFile)

