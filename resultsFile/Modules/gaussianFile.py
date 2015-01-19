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

gaussianFile_defined_vars = [ "date", "version", "machine", "memory", "disk",\
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

class gaussianFile(resultsFile):
   """ Class defining the gaussian file.
   """

   local_vars = list(local_vars)
   defined_vars = list(gaussianFile_defined_vars)

   def get_options(self):
     if self._machine is None:
        self.find_string("\\")
        pos = self._pos
        self.find_next_string("@")
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

   def get_date(self):
     if self._date is None:
       self._date = self.options[0][8]
     return self._date

   def get_author(self):
     if self._author is None:
       self._author = self.options[0][7].lower()
     return self._author

   def get_charge(self):
     if self._charge is None:
       self._charge = float(self.options[3][0].split(',')[0])
     return self._charge

   def get_dipole(self):
     if self._dipole is None:
        self._dipole = []
        dip = self.text[pos]
        for x in dip:
          self._dipole.append(float(x))
     return self._dipole
   
   def get_multiplicity(self):
     if self._multiplicity is None:
       self._multiplicity = int(self.options[3][0].split(',')[1])
     return self._multiplicity

   def get_num_elec(self):
     if self._num_elec is None:
       self._num_elec = self.num_alpha + self.num_beta
     return self._num_elec

   def get_nuclear_energy(self):
     if self._nuclear_energy is None:
       self.find_string("N-N")
       pos = self._pos
       self._nuclear_energy = self.text[pos].replace('=','= ').split()[1]
       self._nuclear_energy = float(self._nuclear_energy.replace('D','E'))
     return self._nuclear_energy

   def get_title(self):
     if self._title is None:
       self._title = self.options[2][0].strip()
     return self._title

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

   def get_cpu_time(self):
     if self._cpu_time is None:
       try:
         self.find_string("Job cpu time")
       except IndexError:
         return None
       pos = self._pos
       line = self.text[pos].split()
       self._cpu_time  = float(line[3])*60*60*24
       self._cpu_time += float(line[5])*60*60
       self._cpu_time += float(line[7])*60
       self._cpu_time += float(line[9])
       self._cpu_time = str(self._cpu_time)+' s'
     return self._cpu_time

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
       end = self._pos
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
          if method.startswith("U"): self._spin_restrict = False
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

   def get_energies(self):
       if self._energies is None:
          self._energies = []
          for m in self.methods:
            if m == 'CASSCF' \
            or m.startswith('R')\
            or m.startswith('U'):
             #if self.point_group == "C01":
                self._energies.append(float(self.options[4][2].split('=')[1]))
             #else:
             #  self._energies.append(float(self.options[4][3].split('=')[1]))
       return self._energies

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

   def get_pot_energies(self):
     if self._pot_energies is None:
        self._pot_energies = []
        for i,e in enumerate(self.Ne_pot_energies):
          self._pot_energies.append (self.ee_pot_energies[i]+\
                          e+self.nuclear_energy)
     return self._pot_energies

   def get_kin_energies(self):
     if self._kin_energies is None:
        self.find_string("N-N")
        pos = self._pos
        self._kin_energies = [float(self.text[pos].replace('=','= ').split()[5])]
     return self._kin_energies

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
              self.find_next_string("X           Y           Z")
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

     try:
         b = self.basis
         for f in b:
           for i in range(len(self._geometry)):
             self._geometry[i].basis = []
             at = self._geometry[i]
             if f.center is at.coord:
               self._geometry[i].basis.append(f)
     except IndexError:
         pass
     return self._geometry

   def get_basis(self):
     if self._basis is None:
        gfprint=False
        gfinput=False
        buffer = map(lambda x: x.lower(),self.options[1][0].split())
        if 'gfprint' in buffer: gfprint=True
        elif 'gfinput' in buffer: gfinput=True

        if gfprint:
          Polar=False
          try:
            self.find_string("Standard basis")
          except:
            try:
              self.find_string("General basis")
            except:
              self.find_string("Basis read from chk")
          pos = self._pos
          if "5D" in self.text[pos]: Polar=True
          try:
            self.find_next_string("AO basis set")
            pos = self._pos
          except IndexError:
            return None
          try:
            self.find_string("Pseudopotential Parameters")
            end = self._pos-1
          except:
            pass
          try:
            self.find_string("Integral buffers")
            end = self._pos
          except:
            pass
          try:
            self.find_string("There are")
            end = self._pos
          except:
            pass
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
                    coef2 = float(line[2].replace('D','E'))
                    bf.append( [expo,coef,coef2] )
                 else:
                    bf.append( [expo,coef] )
                 pos += 1
                 line = self.text[pos].split()
             if len(bf) > 0:
               basis_read.append( [index,sym,bf,iatom] )
             if line[0].startswith('==============='):
               pos = end
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
                mylist = [ "G0", "G1+", "G1-", "G2+", "G2-", "G3+", "G3-",
                           "G4+", "G4-" ]
           elif basis[k][0] == "H":
             if not Polar:
               mylist = [ "ZZZZZ", "YZZZZ", "YYZZZ", "YYYZZ", "YYYYZ", "YYYYY", 
                 "XZZZZ", "XYZZZ", "XYYZZ", "XYYYZ", "XYYYY", "XXZZZ", "XXYZZ", "XXYYZ",
                 "XXYYY", "XXXZZ", "XXXYZ", "XXXYY", "XXXXZ", "XXXXY", "XXXXX" ] 
             else:
                mylist = [ "H0", "H1+", "H1-", "H2+", "H2-", "H3+", "H3-",
                                 "H4+", "H4-", "H5+", "H5-" ]
           elif basis[k][0] == "I":
             if not Polar: 
               mylist = [   "ZZZZZZ",   "YZZZZZ",   "YYZZZZ",   "YYYZZZ",
               "YYYYZZ",   "YYYYYZ",   "YYYYYY",   "XZZZZZ",   "XYZZZZ",   "XYYZZZ",
               "XYYYZZ",   "XYYYYZ",   "XYYYYY",   "XXZZZZ",   "XXYZZZ",   "XXYYZZ",
               "XXYYYZ",   "XXYYYY",   "XXXZZZ",   "XXXYZZ",   "XXXYYZ",   "XXXYYY",
               "XXXXZZ",   "XXXXYZ",   "XXXXYY",   "XXXXXZ",   "XXXXXY",   "XXXXXX" ]
             else:
                mylist = [ "G0", "G1+", "G1-", "G2+", "G2-", "G3+", "G3-",
                                 "H4+", "H4-", "H5+", "H5-", "I6+", "I6-" ]

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
              gauss.sym  = normalize_basis_name(buffer[0])
              if len(c) == 3:
                 if gauss.sym == "s":
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

   def get_mo_sets(self):
      if self._mo_sets is None:
         self._mo_sets = []
         self._mo_types = []
         posend = {}
       # RHF/UHF/ROHF/CASSCF
         method = self.methods[0]
         try:
           if method.startswith("R") \
             or method == 'CASSCF':
             index = method
             posend[index] = [0,0]
             self._mo_types.append(index)
             self.find_string("Molecular Orbital Coefficients")
           elif method.startswith("U"):
             index = 'Alpha'
             posend[index] = [0,0]
             self._mo_types.append(index)
             self.find_string("Alpha Molecular Orbital Coefficients")
           else:
             raise TypeError
           self.find_next_string(" 1         2")
           posend[index] = []
           posend[index].append(self._pos)
           end = self._pos + (4+len(self.basis))/5 * ( len(self.basis) + 3 )
           posend[index].append(end)
           if method.startswith("U"):
             self.find_string("Beta Molecular Orbital Coefficients")
             index = 'Beta'
             self._mo_types.append(index)
             self.find_next_string(" 1         2")
             posend[index] = []
             posend[index].append(self._pos)
             end = self._pos + (4+len(self.basis))/5 * ( len(self.basis) + 3 )
             posend[index].append(end)
         except IndexError:
           pass
         try:
           self.find_string("Natural Orbital Coefficients")
           index = 'Natural'
           self._mo_types.append(index)
           self.find_next_string(" 1         2")
           posend[index] = []
           posend[index].append(self._pos)
           end = self._pos + (4+len(self.basis))/5 * ( len(self.basis) + 2 )
           posend[index].append(end)
         except IndexError:
           pass
         self._mo_sets = {}
         for index  in self._mo_types:
            pos,end = posend[index]
            curatom = 1
            vectors = []
            while pos < end:
               pos += 1
               if index == 'Natural':
                 syms = [ None for k in range(5) ]
               else:
                 line = self.text[pos].split()
                 syms = []
                 for sym in line:
                   if '(' in sym and ')' in sym:
                     syms.append ( sym.split(')')[0].split('(')[1] )
                   else:
                     syms.append ( sym )
                 pos += 1
               line = self.text[pos][19:].replace('-',' -')
               line = line.replace('**********','-99999999.')
               line = line.split()
               if line[0].lower() == "unable":
                 pos +=1 
                 line = line.split()
               eigval = [ float(k) for k in line ]
               pos += 1
             # Build vectors
               begin = pos
               lmax = len(eigval)
               for l in range(lmax):
                 pos = begin
                 v = orbital()
                 v.basis = self.basis
                 v.set = index
                 v.eigenvalue = eigval[l]
                 line = self.text[pos][1:21].split()
                 if len(syms) >0: v.sym = syms[l]
                 while len(line) > 0 and pos < end:
                    line = self.text[pos][21:].strip().replace('-',' -')
                    line = line.split()
                    try:
                      v.vector.append(float(line[l]))
                    except:
                      break
                    pos += 1
                    line = self.text[pos][1:21].strip().split()
                    if "DENSITY" in line:
                      line = []
                 vectors.append(v)
            self._mo_sets[index] = vectors
         if self._mo_sets == {}:
             self._mo_sets = None
      return self._mo_sets
         
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

   def get_num_alpha(self):
      if self._num_alpha is None:
       self.find_string("alpha electrons")
       pos = self._pos
       self._num_alpha = int(self.text[pos].split()[0])
      return self._num_alpha 

   def get_num_beta(self):
      if self._num_beta is None:
       self.find_string("beta electrons")
       pos = self._pos
       self._num_beta = int(self.text[pos].split()[3])
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
      method = self.methods[0]
      if self._csf is None:

        csf = []
        # Mono-determinant case
        if method.startswith("R") and not method.startswith("RO"):
           new_csf = CSF()
           new_spin_det = []
           for i in range(self.num_alpha):
               new_spin_det.append(i)
           new_csf.append(1.,new_spin_det,new_spin_det)
           csf.append(new_csf)
           self._determinants_mo_type = self.mo_types[-1]

        elif method.startswith("RO") or method.startswith("U"):
           new_csf = CSF()
           new_spin_deta = []
           new_spin_detb = []
           for i in range(self.num_alpha):
               new_spin_deta.append(i)
           for i in range(self.num_beta):
               new_spin_detb.append(i)
           new_csf.append(1.,new_spin_deta,new_spin_detb)
           csf.append(new_csf)
           self._determinants_mo_type = self.mo_types[-1]

        # Multi-determinant case
        if method == 'CASSCF':
           try:
             self.find_string('Enter MCSCF program')
             self.find_next_string('CORE-ORBITALS')
             pos = self._pos
             buffer = self.text[pos].split('=')
             ncore = int(buffer[-1])
             self.find_string('no. active ELECTRONS')
             pos = self._pos
             buffer = self.text[pos].split()
             nact  = int (buffer[4])
             self.find_next_string('PRIMARY BASIS FUNCTION')
             pos = self._pos
             mostr = self.text[pos].split('=')[1].split()
             tempcsf_a = []
             tempcsf_b = []
             for i in range(ncore):
               tempcsf_a.append(i)
               tempcsf_b.append(i)
             to_add = tempcsf_a
             old=0
             for i in mostr:
               i = int(i)
               if i <= old: to_add=tempcsf_b
               to_add.append(i+ncore-1)
               old=i
             this_csf = CSF()
             this_csf.append(1.,tempcsf_a,tempcsf_b)
             csf.append(this_csf)
             pos = self._pos+1
             self.find_next_string("NO OF BASIS FUNCTIONS")
             end = self._pos
             while pos < end:
               this_csf = CSF()
               tempcsf_a = []
               tempcsf_b = []
               for tempcsf in [tempcsf_a,tempcsf_b]:
                 pos += 1
                 buffer = self.text[pos].split()
                 for i in range(ncore):
                   tempcsf.append(i)
                 for i in buffer:
                   orb = int (i)+ncore-1
                   tempcsf.append(orb)
               this_csf.append(1.,tempcsf_a,tempcsf_b)
               pos += 1
               csf.append(this_csf)
           except IndexError:
             pass

        if csf != []:
           self._csf = csf
      return self._csf


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

   def get_det_coefficients(self):
      if self._det_coefficients is None:
        # Mono-determinant case
        if self.options[0][4].startswith('R') or\
           self.options[0][4].startswith('U'):
           self._det_coefficients = [ [1.] ]
        # Multi-determinant case TODO
        if self.csf is not None:
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
        # Mono-determinant case
        if self.options[0][4].startswith('R') or\
           self.options[0][4].startswith('U'):
           self._csf_coefficients=[ [1.] ]
        # Multi-determinant case TODO
        elif self.options[0][4] == 'CASSCF':
           self._csf_coefficients = [ [] for i in range(self.num_states) ]
           self._pos = 0
           regexp = re.compile("\( *([0-9]+)\)([ |-][0-9]*\.[0-9]*) *")
           for state in range(self.num_states):
             for i in self.csf:
               self._csf_coefficients[state].append(0.)
             cur_csf=0
             self.find_string("2ND ORD PT ENERGY")
             self._pos += 1
             self.find_next_string("EIGENVALUE ")
             pos = self._pos + 1
             matches = [1]
             while matches != []:
               matches = regexp.findall(self.text[pos])
               for m in matches:
                 self._csf_coefficients[state][int(m[0])-1] = float(m[1])
               pos += 1
      return self._csf_coefficients

   def get_occ_num(self):
     if self._occ_num is None:
      occ = {}
      method = self.methods[0]
      for motype in self.mo_types:
          occ[motype] = [ 0. for mo in self.mo_sets[motype] ]
      if method.startswith("R") and not method.startswith("RO"):
        for i in range(self.num_alpha):
          occ[motype][i] = 2.
      elif method.startswith("RO"):
        for i in range(self.num_beta):
          occ[motype][i] = 2.
        for i in range(self.num_beta,self.num_alpha):
          occ[motype][i] = 1.
      elif self.mulliken_mo is not None:
        for motype in self.mo_types:
            occ[motype] = [ mo.eigenvalue for mo in self.mulliken_mo ]
      self._occ_num = occ
     return self._occ_num 

   def get_num_states(self):
      if self._num_states is None:
        self._num_states=1
        try:
          self.find_string("NUMBER OF STATES REQUESTED")
          pos = self._pos
          buffer = self.text[pos].split('=')
          self._num_states = int(buffer[1])
        except IndexError:
          pass
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


   class two_e_int_array(object):

        def __init__(self,filename):
           self._file = fortranBinary(filename,"rb")

        def get_form(self):
           return self._file.form 

        def set_form(self,form):
           self._file.form = form

        def __iter__(self):
          for rec in self._file:
            Nint = rec.pop(0)
            rec_i = rec[:4*Nint]
            rec_f = rec[4*Nint:]
            del rec
            for m in range(Nint):
              n = 4*m
              i = rec_i[n]
              j = rec_i[n+1]
              k = rec_i[n+2]
              l = rec_i[n+3]
              v = rec_f[m]
              int = integral()
              int.indices = (i, j, k, l)
              int.value = v
              yield int

        def __getitem__(self,index):
          self._file.seek(0)
          rec = self._file[0]
          Nint = rec.pop(0)
          while (index >= Nint):
            index -= Nint
            rec = self._file.next()
            Nint = rec.pop(0)

          rec_i = rec[:4*Nint]
          rec_f = rec[4*Nint:]
          del rec
          m = index
          n = 4*m
          i = rec_i[n]
          j = rec_i[n+1]
          k = rec_i[n+2]
          l = rec_i[n+3]
          v = rec_f[m]
          int = integral()
          int.indices = (i, j, k, l)
          int.value = v
          return int

        form = property ( get_form, set_form )

# Output Routines
# ---------------

def write_vec(res, file):
   for motype in res.mo_types:
     print >>file, " $VEC"
     moset = res.mo_sets[motype]
     for idx,mo in enumerate(moset):
       v = list(mo.vector)
       line = 1
       while (len(v) > 0):
         fields = []
         monum = "%4d"%(idx )
         monum = monum[-2:]
         fields.append ( monum )
         linum = "%4d"%(line)
         linum = linum[-3:] 
         fields.append ( linum )
         for n in range(5):
           fields.append ( "%15.8E"%v.pop(0) )
         print >>file, ''.join(fields)
         line += 1
     print >>file, " $END"

def write_data(res, file):
  print >>file, " $DATA"
  print >>file, res.title
  point_group = res.point_group
  N = point_group[1]
  if N != '1':
     point_group = point_group[0]+'N'+point_group[2:]
     point_group += '    '+N
  print >>file, point_group
  if N != '1':
    print >>file, ""
  for at in res.geometry:
    print >>file, "%10s  %3.1f  %15.8f %15.8f %15.8f" % tuple( \
     [at.name.ljust(10), at.charge]+list(at.coord) )
    for contr in at.basis:
      sym = contr.sym
      if sym == 's' :
         doPrint = True
      elif sym == 'x' :
         sym = 'P'
         doPrint = True
      elif sym == 'xx' :
         sym = 'D'
         doPrint = True
      elif sym == 'xxx' :
         sym = 'F'
         doPrint = True
      elif sym == 'xxxx' :
         sym = 'G'
         doPrint = True
      elif sym == 'xxxxx' :
         sym = 'H'
         doPrint = True
      else:
         doPrint = False
      if doPrint:
         print >>file, "%4s%8d"%(sym, len(contr.prim))
         for i, p in enumerate(contr.prim):
            if not isinstance(p,gaussian):
               raise TypeError("Basis function is not a gaussian")
            print >>file, "%6d%26.10f%12.8f"%(i+1,p.expo,contr.coef[i])
    print >>file, ""
  print >>file, " $END"


      
fileTypes.append(gaussianFile)
  
if __name__ == '__main__':
   main(gaussianFile)
