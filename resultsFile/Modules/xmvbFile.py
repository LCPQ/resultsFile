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

xmvbFile_defined_vars = [ "date", "version", "machine", "memory", "disk",\
                "cpu_time", "author", "title", "units", "methods", "options", \
                "spin_restrict", "conv_threshs", "energies", \
                "ee_pot_energies", \
                "Ne_pot_energies", "pot_energies", \
                "kin_energies", "point_group", "num_elec", \
                "charge", "multiplicity","nuclear_energy","dipole","geometry",\
                "basis","mo_sets","mo_types","mulliken_mo","mulliken_ao",\
                "mulliken_atom","lowdin_ao", "mulliken_atom","lowdin_atom",\
                "two_e_int_ao", "determinants", "num_alpha", "num_beta",\
                "closed_mos", "active_mos", "virtual_mos", \
                "determinants_mo_type", "det_coefficients", \
                "csf_mo_type", "csf_coefficients", "symmetries", "occ_num", \
                "csf", "num_states"]

class xmvbFile(resultsFile):
   """ Class defining the xmvb file.
   """

   local_vars = list(local_vars)
   defined_vars = list(xmvbFile_defined_vars)

   def get_options(self):
     if self._machine is None:
        self.find_string("\\")
        pos = self._pos
        self.find_next_string("\\@")
        end = self._pos
        buffer = ""
        self._options = []
        for line in self.text[pos:end+1]:
           buffer += line[1:].replace('\n','')
        buffer = buffer.split('\\\\')
        for l in buffer:
          self._options.append(l.split('\\'))
        self._options.pop()
     return self._options

   def get_machine(self):
     if self._machine is None:
        self._machine = self.options[0][2][5:].lower()
     return self._machine

   def get_version(self):
     if self._version is None:
       self._version = self.options[4][0].split('=')[1]
     return self._version

   def get_author(self):
     if self._author is None:
       self._author = self.options[0][7].lower()
     return self._author

   def get_charge(self):
     if self._charge is None:
       self._charge = float(self.options[3][0].split(',')[0])
     return self._charge

   def get_disk(self):
     if self._disk is None:
       try:
         self.find_string("File lengths")
       except IndexError:
         return None
       pos = self._pos
       line = self.text[pos].split()
       disk = 0
       for i in line[4::2]:
         disk += float(i)
       disk = disk/1000.
       if disk > 1.:
         self._disk = str(disk)+" Gb"
       else:
         disk *= 1000.
         if disk > 1.:
           self._disk = str(disk)+" Mb"
         else:
           disk *= 1000.
           self._disk = str(disk)+" kb"
     return self._disk

   def get_memory(self):
     if self._memory is None:
       try:
         self.find_string("Leave Link")
       except IndexError:
         return None
       pos = self._pos
       line = self.text[pos].split()
       memory = float(line[10])*8. / 1000000000.
       if memory > 1.:
         self._memory = str(memory)+" Gb"
       else:
         memory *= 1000.
         if memory > 1.:
           self._memory = str(memory)+" Mb"
         else:
           memory *= 1000.
           self._memory = str(memory)+" kb"
     return self._memory

   def get_symmetries(self):
     if self._symmetries is None:
       try:
         self.find_string("There are")
       except IndexError:
         return None
       pos = self._pos
       begin = pos
       try:
         self.find_next_string("Integral")
       except IndexError:
         return None
       end = self._pos-1
       sym = []
       for k in range(begin,end):
         buffer = self.text[k].split()
         sym.append([buffer[8],int(buffer[2])])
       self._symmetries = sym
     return self._symmetries

   def get_units(self):
     if self._units is None:
       try:
         self.find_string("Coordinates")
       except IndexError:
         return None
       pos = self._pos
       units = self.text[pos].split()[4][1:-1]
       if units != 'Angstroms':
          self._units = 'BOHR'
       else:
          self._units = 'ANGS'
     return self._units

   def get_methods(self):
     if self._methods is None:
        methods = []
        methods.append(self.options[0][4])
        self._methods = methods
     return self._methods

   def get_spin_restrict(self):
       if self._spin_restrict is None:
          method = self.methods[0]
          self._spin_restrict = True 
          if method == 'UHF': self._spin_restrict = False
       return self._spin_restrict


   def get_conv_threshs(self):
       if self._conv_threshs is None:
          self._conv_threshs = []
          for m in self.methods:
            if m == 'RHF' or m == 'UHF' or m == 'ROHF':
              self.find_string("SCF Done")
              pos = self._pos + 1
              self._conv_threshs.append(float(self.text[pos].split()[2]))
            if m == 'CASSCF':
              self.find_string("Enter MCSCF program")
              self.find_next_string("USED ACCURACY IN CHECKING CONVEGERGENCE")
              pos = self._pos
              self._conv_threshs.append(float(self.text[pos].split('=')[1]))
          if self._conv_threshs == []:
              self._conv_threshs = None
       return self._conv_threshs

   def get_ee_pot_energies(self):
     if self._ee_pot_energies is None:
        self._ee_pot_energies = [] 
        for i,e in enumerate(self.kin_energies):
          self._ee_pot_energies.append(self.energies[i]\
                          -self.nuclear_energy\
                          -self.kin_energies[i]\
                          -self.Ne_pot_energies[i])
     return self._ee_pot_energies

   def get_Ne_pot_energies(self):
     if self._Ne_pot_energies is None:
        self.find_string("N-N")
        pos = self._pos
        self._Ne_pot_energies = [float(self.text[pos].replace('=','= ').split()[3])]
     return self._Ne_pot_energies

   def get_point_group(self):
     if self._point_group is None:
        self._point_group = self.options[4][-1].split()[0].split('=')[1]
     return self._point_group 

   def get_geometry(self):
     if self._geometry is None:
        self._geometry = []
        self._pos = 0
        pos=0
        try:
           while True:
              pos = self._pos
              self.find_next_string("Number     Number      Type")
              self._pos += 1
        except IndexError:
           pass
        pos +=1
        self._pos=pos
        self.find_next_string("-----")
        end = self._pos
        while pos<end:
           temp = atom()
           buffer = self.text[pos].split()
           temp.charge = float(buffer[1])
           temp.coord  = (float(buffer[3]), float(buffer[4]), float(buffer[5]))
           self._geometry.append(temp)
           pos += 1 
        for i,line in enumerate(self.options[3][1:]):
           buffer = line.split(',')
           self._geometry[i].name = buffer[0]
#       try:
#          b = self.basis
#          for f in b:
#            for at in self._geometry:
#              if f.center is at.coord:
#                at.basis.append(f)
#       except IndexError:
#          pass
     return self._geometry

   def get_basis(self):
     if self._basis is None:
        gfprint=False
        gfinput=False
        buffer = self.options[1][0].split()
        if 'GFPRINT' in buffer: gfprint=True
        elif 'GFINPUT' in buffer: gfinput=True
       
        if gfprint:
          Polar=False
          try:
            self.find_string("Standard basis")
          except:
            self.find_string("General basis")
          pos = self._pos
          if "5D" in self.text[pos]: Polar=True
          try:
            self.find_next_string("AO basis set:")
            pos = self._pos
          except IndexError:
            return None
          self.find_string("Integral buffers")
          end = self._pos
          doLoop=True
          while (doLoop):
            try:
              self.find_prev_string("There are")
              end = self._pos
            except:
              doLoop = False
          pos += 1
          basis_read = []
          line = self.text[pos].split()
          iatom=0
          atom = line[1]
          while pos < end:
             if line[0] == 'Atom':
               index = int(line[3])
               sym = line[4]
               nfunc = int(line[5])
               if atom != line[1]:
                 iatom += 1
               atom = line[1]
               bf = []
               pos+=1
               line = self.text[pos].split()
               for k in range(nfunc):
                 expo = float(line[0].replace('D','E'))
                 coef = float(line[1].replace('D','E'))
                 if sym == "SP":
                    coef2 = float(line[2])
                    bf.append( [expo,coef,coef2] )
                 else:
                    bf.append( [expo,coef] )
                 pos += 1
                 line = self.text[pos].split()
             if len(bf) > 0:
               basis_read.append( [index,sym,bf,iatom] )
        else:
          print "GFPRINT should be present in the gaussian keywords."
          return None
        Nmax = basis_read[len(basis_read)-1][0]
        basis = [None for i in range(Nmax)]
        for b in basis_read:
           basis[b[0]-1] = [b[1],b[2],b[3]]
        NotNone = 0
        ReadNone = False
        for i in range(len(basis)-1,-1,-1):
           if basis[i] == None:
              ReadNone = True
              basis[i] = list(basis[i+NotNone])
           else:
              if ReadNone:
                 NotNone = 0
                 ReadNone = False
              NotNone += 1

        k=0
        while k<len(basis):
           if basis[k][0] == "S":
              mylist = []
           elif basis[k][0] == "P":
              mylist = [ "X", "Y", "Z" ]
           elif basis[k][0] == "SP":
              mylist = [ "S", "X", "Y", "Z" ]
           elif basis[k][0] == "D":
              if not Polar:
                mylist = [ "XX", "YY", "ZZ", "XY", "XZ", "YZ" ]
              else:
                mylist = [ "D0", "D1+", "D1-", "D2+", "D2-" ]
           elif basis[k][0] == "F":
              if not Polar:
                mylist = [ "XXX", "YYY", "ZZZ", "YYX", "XXY", "XXZ", "ZZX", "ZZY",
                                "YYZ", "XYZ" ]
              else:
                mylist = [ "F0", "F1+", "F1-", "F2+", "F2-", "F3+", "F3-" ]
           elif basis[k][0] == "G":
              if not Polar:
                mylist = [ "ZZZZ", "ZZZY", "YYZZ", "YYYZ", "YYYY", "ZZZX",
                                "ZZXY", "YYXZ","YYYX", "XXZZ", "XXYZ", 
                              "XXYY", "XXXZ","XXXY", "XXXX" ]
              else:
                mylist = [ "G0", "G1+", "G1-", "G2+", "G2+", "G3+", "G3-",
                           "G4+", "G4-" ]
          #elif basis[k][0] == "H":
          #   mylist = [ "XXXXX", "YYYYY", "ZZZZZ", "XXXXY", "XXXXZ", "YYYYX", \
          #                   "YYYYZ", "ZZZZX", "ZZZZY", "XXXYY", "XXXZZ", \
          #                   "YYYXX", "YYYZZ", "ZZZXX", "ZZZYY", "XXXYZ", \
          #                   "YYYXZ", "ZZZXY", "XXYYZ", "XXZZY", "YYZZX" ]
          #elif basis[k][0] == "I":
          #   mylist = [ "XXXXXX", "YYYYYY", "ZZZZZZ", "XXXXXY", "XXXXXZ", \
          #                   "YYYYYX", "YYYYYZ", "ZZZZZX", "ZZZZZY", "XXXXYY", \
          #                   "XXXXZZ", "YYYYXX", "YYYYZZ", "ZZZZXX", "ZZZZYY", \
          #                   "XXXXYZ", "YYYYXZ", "ZZZZXY", "XXXYYY", "XXXZZZ", \
          #                   "YYYZZZ", "XXXYYZ", "XXXZZY", "YYYXXZ", "YYYZZX", \
          #                   "ZZZXXY", "ZZZYYX" ]
           mylist = map(normalize_basis_name,mylist)
           for i in mylist[:-1]:
              basis[k][0] = i
              basis.insert(k,list(basis[k]))
              k+=1
           for i in mylist[-1:]:
              basis[k][0] = i
           k+=1

        self._basis = []
        for buffer in basis:
            contr = contraction()
            for c in buffer[1]:
              gauss = gaussian()
              atom = self.geometry[int(buffer[2])]
              gauss.center = atom.coord
              gauss.expo = c[0]
              gauss.sym  = buffer[0]
              if len(c) == 3:
                 if gauss.sym == "S":
                    contr.append(c[1],gauss)
                 else:
                    contr.append(c[2],gauss)
              else:
                 contr.append(c[1],gauss)
            self._basis.append(contr)

     return self._basis

   def get_mo_types(self):
      if self._mo_types is None:
        self.get_mo_sets()
      return self._mo_types 

   def get_mulliken_mo(self):
      if self._mulliken_mo is None:
         pass
      return self._mulliken_mo

   def get_mulliken_ao(self):
       if self._mulliken_ao is None:
          pass
       return self._mulliken_ao

   def get_lowdin_ao(self):
       if self._lowdin_ao is None:
          pass
       return self._lowdin_ao


   def get_mulliken_atom(self):
       if self._mulliken_atom is None:
          try:
             self.find_string("Mulliken atomic charges")
          except IndexError:
             return None
          self._pos += 2
          pos = self._pos
          self.find_next_string("Sum of Mulliken")
          end = self._pos
          line = self.text[pos].split()
          vm = orbital()
          vm.set = "Mulliken_Atom"
          vm.eigenvalue = 0.
          while pos < end:
             value = float(line[2])
             vm.vector.append(value)
             vm.eigenvalue += value
             pos += 1
             line = self.text[pos].split()
          self._mulliken_atom = vm
       return self._mulliken_atom

   def get_lowdin_atom(self):
       if self._lowdin_atom is None:
          pass
       return self._lowdin_atom

   def get_two_e_int_ao(self):
      if self._two_e_int_ao is None:
         pass
      return self._two_e_int_ao

   def get_closed_mos(self):
      if self._closed_mos is None:
         result = []
         maxmo = len(self.mo_sets[self.determinants_mo_type])
         for orb in range(maxmo):
           present = True
           for det in self.determinants:
               for spin_det in det:
                   if orb not in det[spin_det]: present = False
           if present: result.append(orb)
           self._closed_mos = result
      return self._closed_mos

   def get_virtual_mos(self):
      if self._virtual_mos is None:
         result = []
         minmo = len(self.closed_mos)
         maxmo = len(self.mo_sets[self.determinants_mo_type])
         for orb in range(minmo,maxmo):
           present = False
           for det in self.determinants:
               for spin_det in det:
                   if orb in det[spin_det]: present = True
           if not present: result.append(orb)
           self._virtual_mos = result
      return self._virtual_mos

   def get_active_mos(self):
      if self._active_mos is None:
         cl = self.closed_mos
         vi = self.virtual_mos
         maxmo = len(self.mo_sets[self.determinants_mo_type])
         result = []
         for i in range(maxmo):
           present = i in cl or i in vi
           if not present:
             result.append(i)
         self._active_mos = result
      return self._active_mos

   def get_occ_num(self):
     if self._occ_num is None:
      if self.mulliken_mo is not None:
        occ = {}
        for motype in self.mo_types:
            occ[motype] = [ mo.eigenvalue for mo in self.mulliken_mo ]
        if occ != {}:
          self._occ_num = occ
     return self._occ_num 

   def get_date(self):
     if self._date is None:
       self.find_string("Job started at")
       pos = self._pos
       self._date = self.text[pos][16:]
     return self._date

   def get_multiplicity(self):
     if self._multiplicity is None:
       self.find_string("nmul=")
       pos = self._pos
       self._multiplicity = int(self.text[pos].split('=')[-1])
     return self._multiplicity

   def get_num_elec(self):
     if self._num_elec is None:
       self.find_string("nelectron=")
       pos = self._pos
       self._num_elec= int(self.text[pos].split('=')[1].split()[0])
     return self._num_elec

   def get_nuclear_energy(self):
     if self._nuclear_energy is None:
       self.find_string("Nuclear Repulsion Energy:")
       pos = self._pos
       self._nuclear_energy = float(self.text[pos].split(':')[1])
     return self._nuclear_energy

   def get_title(self):
     if self._title is None:
       self.find_string("End of Input")
       pos = self._pos+3
       self._title = self.text[pos].strip()
     return self._title

   def get_cpu_time(self):
     if self._cpu_time is None:
       try:
         self.find_string("Cpu for the Job:")
       except IndexError:
         return None
       pos = self._pos
       self._cpu_time = self.text[pos].split(':')[1].split()[0]+' s'
     return self._cpu_time

   def get_dipole(self):
     if self._dipole is None:
        self.find_string("Dipole moment")
        pos = self._pos+2
        line = self.text[pos].split()
        self._dipole = []
        self._dipole.append(float(line[1]))
        self._dipole.append(float(line[3]))
        self._dipole.append(float(line[5]))
     return self._dipole


   def get_energies(self):
       if self._energies is None:
          self._energies = []
          self.find_string("Total Energy:")
          pos = self._pos
          self._energies.append(float(self.text[pos].split(':')[1]))
       return self._energies

   def get_pot_energies(self):
     if self._pot_energies is None:
        self._pot_energies = [] 
        self.find_string("Potential energy:")
        pos = self._pos
        self._pot_energies.append(float(self.text[pos].split(':')[1]))
     return self._pot_energies

   def get_kin_energies(self):
     if self._kin_energies is None:
        self._kin_energies = [] 
        self.find_string("Kinetic energy:")
        pos = self._pos
        self._kin_energies.append(float(self.text[pos].split(':')[1]))
     return self._kin_energies

   def get_mo_sets(self):
      if self._mo_sets is None:
         self._mo_types = ['VB']
         posend = {}
         self.find_string("norb=")
         norb = int(self.text[self._pos].split("norb=")[1].split()[0])
         self.find_string("OPTIMIZED ORBITALS")
         pos = self._pos+3
         iorb=0
         vectors = []
         while(iorb<norb):
            lorb = len(self.text[pos].split())
            pos += 1
            begin = pos
            for l in range(lorb):
              iorb+=1
              pos = begin
              line = self.text[pos].split()
              v = orbital()
              v.set = self.mo_types[0]
              v.basis = None
              v.eigenvalue = float(iorb)
              while len(line) > 0:
                v.vector.append(float(line[l+1]))
                pos += 1
                line = self.text[pos].split()
              vectors.append(v)
         self._mo_sets = {}
         self._mo_sets['VB'] = vectors
      return self._mo_sets
         
   def get_num_alpha(self):
      if self._num_alpha is None:
         self._num_alpha = self.num_elec/2 + (self.multiplicity-1)/2
      return self._num_alpha

   def get_num_beta(self):
      if self._num_beta is None:
         self._num_beta = self.num_elec/2 - (self.multiplicity-1)/2
      return self._num_beta

   def get_determinants_mo_type(self):
      if self._determinants_mo_type is None:
         self._determinants_mo_type = self.mo_types[-1]
      return self._determinants_mo_type

   def get_csf_mo_type(self):
      if self._csf_mo_type is None:
         self._csf_mo_type = self.determinants_mo_type
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
      if self._csf is None:
        csf = []
        csf_coefficients = []
        self.find_string('COEFFICIENTS OF DETERMINANTS')
        self.find_next_string('1')
        pos = self._pos
        buffer = self.text[pos]
        while buffer.strip() != "":
          tempcsf_a = []
          tempcsf_b = []
          buffer = self.text[pos]
          coef = float(buffer.split()[1])
          ii=0
          while ii < self.num_elec:
            buffer = buffer[24:].split()
            for i in buffer:
              ii+=1
              if ii <= self.num_alpha:
                tempcsf_a.append(int(i)-1)
              else:
                tempcsf_b.append(int(i)-1)
            pos += 1
            buffer = self.text[pos]
          this_csf = CSF()
          this_csf.append(1.,tempcsf_a,tempcsf_b)
          csf.append(this_csf)
          csf_coefficients.append([coef])
        if csf != []:
           self._csf = csf
           self._csf_coefficients = csf_coefficients
      return self._csf


   def get_det_coefficients(self):
      if self._det_coefficients is None:
        if self.csf is not None:
           self._det_coefficients = []
           csf = self.csf
           vector = []
           for state_coef in self.csf_coefficients:
              for i,c in enumerate(state_coef):
                 for d in csf[i].coefficients:
                    vector.append(c*d)
              self._det_coefficients.append(vector)
      return self._det_coefficients

   def get_csf_coefficients(self):
      if self._csf_coefficients is None:
        self.get_csf()
      return self._csf_coefficients

   def get_num_states(self):
      if self._num_states is None:
        self._num_states=1
      return self._num_states

   to_remove = []
   for i, j in local_vars:
     if i in resultsFile_defined_vars:
        to_remove.append( (i,j) )
   for i in to_remove:
      local_vars.remove(i)

   for i, j in local_vars:
     if i not in defined_vars:
        exec build_get_funcs(i) in locals()
     exec build_property(i,j) in locals()
   del to_remove, i, j


fileTypes.insert(0,xmvbFile)
  
if __name__ == '__main__':
   main(xmvbFile)

###### END #####

