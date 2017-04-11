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

gamessFile_defined_vars = [ "date", "version", "machine", "memory", "disk",\
                "cpu_time", "author", "title", "units", "methods", "options", \
                "spin_restrict", "conv_threshs", "energies", \
                "one_e_energies", "two_e_energies", "ee_pot_energies", \
                "Ne_pot_energies", "pot_energies", \
                "kin_energies", "virials", "point_group", "num_elec", \
                "charge", "multiplicity","nuclear_energy","dipole","geometry",\
                "basis","mo_sets","mo_types","mulliken_mo","mulliken_ao",\
                "mulliken_atom","lowdin_ao", "mulliken_atom","lowdin_atom",\
                "two_e_int_ao", "determinants", "num_alpha", "num_beta",\
                "closed_mos", "active_mos", "virtual_mos", \
                "determinants_mo_type", "det_coefficients", \
                "csf_mo_type", "csf_coefficients", "symmetries", "occ_num", \
                "csf", "num_states", "two_e_int_ao_filename", "pseudo", 
                "one_e_int_ao_filename", "atom_to_ao_range", "gradient_energy" ]

class gamessFile(resultsFile):
   """ Class defining the gamess file.
   """

   local_vars = list(local_vars)
   defined_vars = list(gamessFile_defined_vars)

   exec get_data('date',"EXECUTION OF GAMESS BEGUN",'4:') in locals()
   exec get_data('machine',"Files used on the master node",'6:7') in locals()
   exec get_data('version',"GAMESS VERSION",'4:8') in locals()
   exec get_data('num_proc',"MEMDDI DISTRIBUTED OVER",'3:4',"int") in locals()
   exec get_data('num_elec',"NUMBER OF ELECTRONS", '4:5',"int") in locals()
   exec get_data('charge',"CHARGE OF MOLECULE", '4:5',"float") in locals()
   exec get_data('nuclear_energy',"THE NUCLEAR REPULSION ENERGY", '5:6',"float") in locals()

   def get_multiplicity(self):
     if self._multiplicity is None:
       try:
         self.find_string("SPIN MULTIPLICITY")
       except IndexError:
         try:
           self.find_string("STATE MULTIPLICITY")
         except IndexError:
           return None
       pos = self._pos
       if pos is not None:
          line = self.text[pos].split()
          self._multiplicity = int(' '.join(line[3:4]))
     return self._multiplicity

   def get_two_e_int_mo_filename(self):
     try:
       getattr(self,'_two_e_int_mo_filename')
     except AttributeError:
       self._two_e_int_mo_filename = None
     if self._two_e_int_mo_filename is None:
       filename = self.filename.split('.')
       filename.pop()
       filename.append("F08")
       filename = '.'.join(filename)
       self._two_e_int_mo_filename = filename
     return self._two_e_int_mo_filename

   def set_two_e_int_mo_filename(self,value):
     self._two_e_int_mo_filename = value

   def get_two_e_int_ao_filename(self):
     try:
       getattr(self,'_two_e_int_ao_filename')
     except AttributeError:
       self._two_e_int_ao_filename = None
     if self._two_e_int_ao_filename is None:
       filename = self.filename.split('.')
       filename.pop()
       filename.append("F08")
       filename = '.'.join(filename)
       self._two_e_int_ao_filename = filename
     return self._two_e_int_ao_filename

   def set_two_e_int_ao_filename(self,value):
     self._two_e_int_ao_filename = value

   def get_disk(self):
     if self._disk is None:
       try:
         self.find_string("Files used on the master node")
       except IndexError:
         return None
       pos = self._pos
       pos += 1
       line = self.text[pos].split()
       disk = 0
       while len(line) == 8:
         disk += float(line[4])
         pos += 1
         line = self.text[pos].split()
       disk = disk/1000000000.
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
         self.find_string("EXECUTION OF GAMESS TERMINATED NORMALLY")
       except IndexError:
         return None
       pos = self._pos
       pos -= 3
       line = self.text[pos].split()
       self._cpu_time = line[-4]+" s"
     return self._cpu_time

   def get_memory(self):
     if self._memory is None:
       try:
         self.find_string("EXECUTION OF GAMESS TERMINATED NORMALLY")
       except IndexError:
         return None
       pos = self._pos
       pos -= 1
       line = self.text[pos].split()
       memory = float(line[0])*8. / 1000000000.
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

   def get_author(self):
     if self._author is None:
       try:
         self.find_string("Files used on the master node")
       except IndexError:
         return None
       pos = self._pos
       pos += 1
       line = self.text[pos].split()
       self._author = line[2]
     return self._author

   def get_title(self):
     if self._title is None:
       try:
         self.find_string("RUN TITLE")
       except IndexError:
         return None
       pos = self._pos
       pos += 2
       self._title = self.text[pos].strip()
     return self._title

   def get_symmetries(self):
     if self._symmetries is None:
       try:
         self.find_string("DIMENSIONS OF THE SYMMETRY SUBSPACES")
       except IndexError:
         return None
       pos = self._pos
       begin = pos + 1
       try:
         self.find_next_string("DONE")
       except IndexError:
         return None
       end = self._pos-1
       sym = []
       for k in range(begin,end):
         buffer = self.text[k].replace('=',' ').split()
         for i,j in zip(buffer[0::2],buffer[1::2]):
           sym.append([i,int(j)])
       self._symmetries = sym
     return self._symmetries

   def get_units(self):
     if self._units is None:
       try:
         self.find_string("ATOM      ATOMIC                      COORDINATES")
       except IndexError:
         return None
       pos = self._pos
       line = list(self.text[pos].split().pop())
       line.pop()
       line.pop(0)
       self._units = ''.join(line)
     return self._units

   def get_methods(self):
     if self._methods is None:
        methods = []
        options = self.options
        for m in ['SCFTYP', 'CITYP']:
          if options[m] != 'NONE':
              methods.append(options[m])
        self._methods = methods
     return self._methods

   def get_occ_num(self):
     if self._occ_num is None:
        occ = {}
        for motype in self.mo_types:
            occ[motype] = [ 0. for mo in self.mo_sets[motype] ]
            if motype == "Natural":
              for i,mo in enumerate(self.mo_sets[motype]):
                  occ[motype][i] = mo.eigenvalue
            else:
              if self.mulliken_mo is not None:
                for i,mo in enumerate(self.mulliken_mo):
                  occ[motype][i] = mo.eigenvalue
        if occ != {}:
          self._occ_num = occ
     return self._occ_num 

   def get_options(self):
     if self._options is None:
        options = {}
        list_of_keys = ["$CONTRL OPTIONS","INTEGRAL TRANSFORMATION OPTIONS", \
                        "INTEGRAL INPUT OPTIONS","GUESS OPTIONS"]
        for key in list_of_keys:
          try:
            self.find_string(key)
            pos = self._pos + 2
            line = self.text[pos].replace('=',' ')
            tokens =  line.split()[0::2]
            values =  line.split()[1::2]
            while len(tokens) > 0:
              for i,t in enumerate(tokens):
                options[t] = values[i]
              pos += 1
              line = self.text[pos].replace('FRIEND=      ','FRIEND=NONE')
              line = line.replace('=',' ')
              tokens =  line.split()[0::2]
              values =  line.split()[1::2]
              if len(tokens) != len(values):
                print "error: ",line
                sys.exit(0)
          except IndexError:
            pass
        try:
           self.find_string("$SYSTEM OPTIONS")
           pos = self._pos + 2
           line = self.text[pos].split()
           options['MEMORY'] = line[2]
           pos += 1
           line = self.text[pos].split()
           options['MEMDDI'] = line[2]
           pos += 3
           line = self.text[pos].replace('=',' ').split()
           options[line[0]] = line[1]
           pos += 1
           line = self.text[pos].replace('=',' ')
           tokens =  line.split()[0::2]
           values =  line.split()[1::2]
           for i,t in enumerate(tokens):
             options[t] = values[i]
        except IndexError:
           pass
        try:
           self.find_string("PROPERTIES INPUT")
           pos = self._pos + 4
           line = self.text[pos].replace('=',' ')
           tokens =  line.split()[0::2]
           values =  line.split()[1::2]
           while tokens[0] != "EXTRAPOLATION":
             for i,t in enumerate(tokens):
               options[t] = values[i]
             pos += 1
             line = self.text[pos].replace('=',' ')
             tokens =  line.split()[0::2]
             values =  line.split()[1::2]
        except IndexError:
           pass
        try:
           self.find_string(" SCF CALCULATION")
           pos = self._pos + 4
           line = self.text[pos].replace('=',' ')
           tokens =  line.split()[0::2]
           values =  line.split()[1::2]
           while tokens[0] != "DENSITY":
             for i,t in enumerate(tokens):
               options[t] = values[i]
             pos += 1
             line = self.text[pos].replace('=',' ')
             tokens =  line.split()[0::2]
             values =  line.split()[1::2]
           line = self.text[pos].replace('=',' ').split()
           options[line[2]] = line[3]
           pos += 1
           line = self.text[pos].replace('=',' ').split()
           options[line[6]] = line[7]
        except IndexError:
           pass
        try:
           self.find_string("GUGA DISTINCT ROW TABLE")
           pos = self._pos + 3
           line = self.text[pos].replace('=',' ')
           tokens =  line.split()[0::2]
           values =  line.split()[1::2]
           while tokens[0] != "-CORE-":
             for i,t in enumerate(tokens):
               options[t] = values[i]
             pos += 1
             line = self.text[pos].replace('=',' ')
             tokens =  line.split()[0::2]
             values =  line.split()[1::2]
        except IndexError:
           pass
        try:
           self.find_string("GUGA DISTINCT ROW TABLE")
           self.find_next_string("-CORE-")
           pos = self._pos + 1
           line = self.text[pos].replace('=',' ')
           tokens =  line.split()[0::2]
           values =  line.split()[1::2]
           while tokens[0] != "THE":
             for i,t in enumerate(tokens):
               options[t] = values[i]
             pos += 1
             line = self.text[pos].replace('=',' ')
             tokens =  line.split()[0::2]
             values =  line.split()[1::2]
        except IndexError:
           pass
        if options != {}:
          self._options = options
     return self._options

   def get_spin_restrict(self):
       if self._spin_restrict is None:
          options = self.options
          if options['SCFTYP'] == 'RHF':
             self._spin_restrict = True
          elif options['SCFTYP'] == 'ROHF':
             self._spin_restrict = True
          elif options['SCFTYP'] == 'UHF':
             self._spin_restrict = False
          elif options['CITYP'] != 'NONE' or options['SCFTYP'] == 'MCSCF':
             self._spin_restrict = True
       return self._spin_restrict


   def get_conv_threshs(self):
       if self._conv_threshs is None:
          options = self.options
          self._conv_threshs = []
          for m in self.methods:
            if m == 'RHF' or m == 'UHF' or m == 'ROHF':
              self._conv_threshs.append(float(options['CONV']))
            if m == 'GUGA':
              try:
                self.find_string("DAVIDSON METHOD")
                self.find_next_string("CONVERGENCE CRITERION")
                res = float(self.text[self._pos].split('=')[1])
                self._conv_threshs.append(res)
              except IndexError:
                pass
       return self._conv_threshs

   def get_gradient_energy(self):
     if self._gradient_energy is None:
       try:
         self.find_string("GRADIENT OF THE ENERGY")
         self.find_next_string("E'X")
         pos = self._pos+1
         self._gradient_energy = []
         for i in self.geometry:
           buffer = self.text[pos].split()
           label = buffer[1]
           assert label == i.name
           f = map(float,buffer[2:])
           self._gradient_energy.append(f)
           pos +=1 
       except IndexError:
         self._gradient_energy = []
     return self._gradient_energy


   def get_energies(self):
       if self._energies is None:
          self._energies = []
          self._one_e_energies = []
          self._two_e_energies = []
          self._ee_pot_energies = []
          self._Ne_pot_energies = []
          self._pot_energies = []
          self._kin_energies = []
          self._virials = []
          poslist = []
          for m in self.methods:
            if m == 'RHF' or m == 'ROHF':
              try:
                self.find_string("FINAL R")
                poslist.append(self._pos)
              except IndexError:
                pass
            elif m == 'UHF':
              try:
                self.find_string("FINAL U")
                poslist.append(self._pos)
              except IndexError:
                pass
            elif m == 'GUGA':
              try:
                self.find_string("END OF DENSITY MATRIX CALCULATION")
                poslist.append(self._pos)
              except IndexError:
                pass
            elif m == 'MCSCF':
              try:
                self.find_string("FINAL MCSCF ENERGY IS")
                poslist.append(self._pos)
              except IndexError:
                pass
            else:
              try:
                self.find_string("FINAL ENERGY IS")
                poslist.append(self._pos)
              except IndexError:
                pass
            for pos in poslist:
              try:
                self._pos = pos
                self.find_next_string("ENERGY COMPONENTS")
                self.find_next_string("ONE ELECTRON ENERGY")
                pos = self._pos
                line = self.text[pos].split()
                self._one_e_energies.append(float(line[4]))
                self.find_next_string("TWO ELECTRON ENERGY")
                pos = self._pos
                line = self.text[pos].split()
                self._two_e_energies.append(float(line[4]))
                self.find_next_string("TOTAL ENERGY")
                pos = self._pos
                line = self.text[pos].split()
                self._energies.append(float(line[3]))
                self.find_next_string("ELECTRON-ELECTRON POT")
                pos = self._pos
                line = self.text[pos].split()
                self._ee_pot_energies.append(float(line[4]))
                self.find_next_string("NUCLEUS-ELECTRON POT")
                pos = self._pos
                line = self.text[pos].split()
                self._Ne_pot_energies.append(float(line[4]))
                self.find_next_string("NUCLEUS-NUCLEUS POT")
                pos = self._pos
                line = self.text[pos].split()
                self._pot_energies.append(float(line[4]))
                self.find_next_string("TOTAL KINETIC ENE")
                pos = self._pos
                line = self.text[pos].split()
                self._kin_energies.append(float(line[4]))
                self.find_next_string("VIRIAL RATIO")
                pos = self._pos
                line = self.text[pos].split()
                self._virials.append(float(line[4]))
              except IndexError:
                pass
       return self._energies

   def get_one_e_energies(self):
     if self._one_e_energies is None:
        self.get_energies()
     return self._one_e_energies

   def get_two_e_energies(self):
     if self._two_e_energies is None:
        self.get_energies()
     return self._two_e_energies

   def get_ee_pot_energies(self):
     if self._ee_pot_energies is None:
        self.get_energies()
     return self._ee_pot_energies

   def get_Ne_pot_energies(self):
     if self._Ne_pot_energies is None:
        self.get_energies()
     return self._Ne_pot_energies

   def get_pot_energies(self):
     if self._pot_energies is None:
        self.get_energies()
     return self._pot_energies

   def get_kin_energies(self):
     if self._kin_energies is None:
        self.get_energies()
     return self._kin_energies

   def get_virials(self):
     if self._virials is None:
        self.get_energies()
     return self._virials

   def get_point_group(self):
     if self._point_group is None:
        try:
           self.find_string("THE POINT GROUP IS")
           pos = self._pos
           line = self.text[pos].replace(',',' ').split()
           group = line[4]
           self._point_group = group.replace('N',line[6])
        except IndexError:
           pass
     return self._point_group 

   def get_geometry(self):
     if self._geometry is None:
        try:
           self.find_string(" ATOM      ATOMIC                      COORDINATES")
           pos = self._pos
           pos += 2
           self._geometry = []
           buffer = self.text[pos].split()
           while len(buffer) > 1:
              temp = atom()
              temp.name   = buffer[0]
              temp.charge = float(buffer[1])
              temp.coord  = (float(buffer[2]), float(buffer[3]), float(buffer[4]))
              temp.basis  = []
              self._geometry.append(temp)
              pos += 1
              buffer = self.text[pos].split()
        except IndexError:
           pass
        b = self.basis
        if b is not None:
           try:
             for f in b:
               for at in self._geometry:
                 if f.center is at.coord:
                   at.basis.append(f)
           except IndexError:
              pass
     return self._geometry

   def get_dipole(self):
     if self._dipole is None:
        try:
           self.find_string("ELECTROSTATIC MOMENTS")
           pos = self._pos
           pos += 6
           self._dipole = []
           for d in self.text[pos].split():
             self._dipole.append( float(d) )
        except IndexError:
           pass
     return self._dipole
   
   def get_basis(self):
     if self._basis is None:
        try:
           self.find_string("SHELL TYPE")
        except IndexError:
           return None
        pos = self._pos
        try:
           self.find_next_string("TOTAL NUMBER OF")
        except IndexError:
           return None
        end = self._pos
        pos += 4
        basis_read = []
        while pos < end:
           line = self.text[pos].split()
           bf = []
           while len(line) > 1:
              index = int(line[0])
              sym = line[1]
              expo = float(line[3])
              coef = float(line[4])
              if sym == "L":
                 coef2 = float(line[5])
                 bf.append( [expo,coef,coef2] )
              else:
                 bf.append( [expo,coef] )
              pos += 1
              line = self.text[pos].split()
           if len(bf) > 0:
              basis_read.append( [index, sym, bf] )
           pos += 1

        Nmax = basis_read[len(basis_read)-1][0]
        basis = [None for i in range(Nmax)]
        for b in basis_read:
           basis[b[0]-1] = [b[1],b[2]]
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
              mylist = ['s']
           elif basis[k][0] == "P":
              mylist = [ "x", "y", "z" ]
           elif basis[k][0] == "L":
              mylist = [ "s", "x", "y", "z" ]
           elif basis[k][0] == "D":
              mylist = [ "xx", "yy", "zz", "xy", "xz", "yz" ]
           elif basis[k][0] == "F":
              mylist = [ "xxx", "yyy", "zzz", "xxy", "xxz", "xyy", "yyz", "xzz", \
                              "yzz", "xyz" ]
           elif basis[k][0] == "G":
              mylist = [ "xxxx", "yyyy", "zzzz", "xxxy", "xxxz", "xyyy", "yyyz", \
                              "xzzz", "yzzz", "xxyy", "xxzz", "yyzz", "xxyz", \
                              "xyyz", "xyzz" ]
           elif basis[k][0] == "H":
              mylist = [ "xxxxx", "yyyyy", "zzzzz", "xxxxy", "xxxxz", "xyyyy", \
                              "yyyyz", "zzzzx", "yzzzz", "xxxyy", "xxxzz", \
                              "xxyyy", "yyyzz", "xxzzz", "yyzzz", "xxxyz", \
                              "xyyyz", "xyzzz", "xxyyz", "xxyzz", "xyyzz" ]
           elif basis[k][0] == "I":
              mylist = [ "xxxxxx", "yyyyyy", "zzzzzz", "xxxxxy", "xxxxxz", \
                              "xyyyyy", "yyyyyz", "xzzzzz", "zzzzzy", "xxxxyy", \
                              "xxxxzz", "xxyyyy", "yyyyzz", "xxzzzz", "yyzzzz", \
                              "xxxxyz", "xyyyyz", "xyzzzz", "xxxyyy", "xxxzzz", \
                              "yyyzzz", "xxxyyz", "xxxyzz", "xxyyyz", "xyyyzz", \
                              "xxyzzz", "xyyzzz" ]
           map(normalize_basis_name,mylist)
           for i in mylist[:-1]:
              basis[k][0] = i
              basis.insert(k,list(basis[k]))
              k+=1
           for i in mylist[-1:]:
              basis[k][0] = i
           k+=1
        for a,i in enumerate(self.atom_to_ao_range):
           for j in range(i[0]-1,i[1]):
              basis[j].append(a)

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
                 if gauss.sym == "s":
                    contr.append(c[1],gauss)
                 else:
                    contr.append(c[2],gauss)
              else:
                 contr.append(c[1],gauss)
#           contr.normalize()
            self._basis.append(contr)

     return self._basis

   def get_mo_types(self):
      if self._mo_types is None:
        self.get_mo_sets()
      return self._mo_types 

   def get_atom_to_ao_range(self):
      try:
       getattr(self,'_atom_to_ao_range')
      except AttributeError:
       self._atom_to_ao_range = None
      if self._atom_to_ao_range is None:
         self._atom_to_ao_range = []
         posend = []
         try:
             index = self.options['SCFTYP']
             self.find_string(" CONVERGED")
             self.find_next_string("EIGENVECTORS")
             self.find_next_string("1          2")
             posend = []
             posend.append(self._pos)
             self.find_next_string("END OF ")
             posend.append(int(self._pos))
         except IndexError:
             pass
         # GUGA
         try:
             if self.options['CITYP'] == 'NONE' or self.options['SCFTYP'] == 'MCSCF':
                raise IndexError
             self.find_string("INITIAL GUESS ORBITALS")
             self.find_next_string("1          2")
             posend = []
             posend.append(self._pos)
             self.find_next_string("END OF INITIAL ORBITAL SELECTION")
             posend.append(self._pos)
         except IndexError:
             pass
         try:
             self.find_string("NATURAL ORBITALS IN ATOMIC ORBITAL BASIS")
             self.find_next_string("1          2")
             posend = []
             posend.append(self._pos)
             self.find_next_string("END OF DENSITY MATRIX CALCULATION")
             posend.append(self._pos)
         except IndexError:
             pass
         # MCSCF
       # try:
       #     if self.options['SCFTYP'] != 'MCSCF':
       #        raise IndexError
       #     self.find_string("MCSCF NATURAL ORBITALS")
       #     self.find_next_string("1          2")
       #     posend = []
       #     posend.append(self._pos)
       #     self.find_next_string("LZ VALUE ANALYSIS")
       #     posend.append(self._pos)
       # except IndexError:
       #     pass
         try:
             if self.options['SCFTYP'] != 'MCSCF':
                raise IndexError
             self.find_string("MCSCF OPTIMIZED ORBITALS")
             self.find_next_string("1          2")
             posend = []
             posend.append(self._pos)
             self.find_next_string(".....DONE WITH MCSCF ITERATIONS....")
             posend.append(self._pos)
         except IndexError:
             pass
         pos,end = posend
         pos -= 1
         curatom = 1
         pos += 2
         line = self.text[pos].split()
         if len(line) == 0:
            pos += 1
            line = self.text[pos].split()
         pos += 2
         line = self.text[pos].split()
         ao_range = [1,0]
         # Build Atom to AO range
         shift = 0
         is_zero = False
         while len(line) > 0 and pos < end:
               try:
                  il2 = int(line[2])+shift
               except ValueError:
                  line.insert(2,str(curatom))
                  self.text[pos] = ' '.join(line)
                  il2 = int(line[2])+shift
               if il2 == 1:
                 is_zero = False
               elif not is_zero and il2 == 0:
                 il2 -= shift
                 shift += 100
                 il2 += shift
                 is_zero = True
               if curatom < il2:
                  self._atom_to_ao_range.append(list(ao_range))
                  ao_range[1] += 1
                  curatom = il2
                  ao_range[0] = int(line[0])
               else:
                  ao_range[1] = int(line[0])
               pos += 1
               line = self.text[pos].split()
         self._atom_to_ao_range.append(ao_range)
         self._atom_to_ao_range = list(self._atom_to_ao_range[:len(self.geometry)])
      return self._atom_to_ao_range


   def get_mo_sets(self):
      if self._mo_sets is None:
         self._mo_sets = []
         self._mo_types = []
         posend = {}
       # UHF
         if not self.spin_restrict:
           try:
             index = self.options['SCFTYP']+'a'
             self.find_string(" CONVERGED")
             self.find_next_string("EIGENVECTORS")
             self.find_next_string("1          2")
             posend[index] = []
             posend[index].append(self._pos)
             self.find_next_string("-- BETA SET --")
             posend[index].append(int(self._pos)-1)
             self._mo_types.append(index)
             index = self.options['SCFTYP']+'b'
             self.find_next_string("EIGENVECTORS")
             self.find_next_string("1          2")
             posend[index] = []
             posend[index].append(self._pos)
             self.find_next_string("END OF "+index[:3]+" CALCULATION")
             posend[index].append(int(self._pos))
             self._mo_types.append(index)
           except IndexError:
             pass
       # RHF/ROHF
         else:
           try:
             index = self.options['SCFTYP']
             self.find_string(" CONVERGED")
             self.find_next_string("EIGENVECTORS")
             self.find_next_string("1          2")
             posend[index] = []
             posend[index].append(self._pos)
             self.find_next_string("END OF ")
             posend[index].append(int(self._pos))
             self._mo_types.append(index)
           except IndexError:
             pass
           try:
             index = self.options['LOCAL']
             self.find_string("THE BOYS LOCALIZED ORBITALS ARE")
             self.find_next_string("1          2")
             posend[index] = []
             posend[index].append(self._pos)
             self.find_next_string("END OF ")
             posend[index].append(int(self._pos))
             self._mo_types.append(index)
           except IndexError:
             pass
         # GUGA
           try:
             if self.options['CITYP'] == 'NONE' and self.options['SCFTYP'] != 'MCSCF':
                raise IndexError
             index = self.options['CITYP']
             self.find_string("INITIAL GUESS ORBITALS")
             self.find_next_string("1          2")
             posend[index] = []
             posend[index].append(self._pos)
             self.find_next_string("END OF INITIAL ORBITAL SELECTION")
             posend[index].append(self._pos)
             self._mo_types.append(index)
           except IndexError:
             pass
           try:
             index = 'Natural'
             self.find_string("NATURAL ORBITALS IN ATOMIC ORBITAL BASIS")
             self.find_next_string("1          2")
             posend[index] = []
             posend[index].append(self._pos)
             self.find_next_string("END OF DENSITY MATRIX CALCULATION")
             posend[index].append(self._pos)
             self._mo_types.append(index)
           except IndexError:
             pass
         # MCSCF
           try:
             if self.options['SCFTYP'] != 'MCSCF':
                raise IndexError
             index = self.options['SCFTYP']
             self.find_string("MCSCF OPTIMIZED ORBITALS")
             self.find_next_string("1          2")
             posend[index] = []
             posend[index].append(self._pos)
             self.find_next_string("DONE WITH MCSCF ITERATIONS")
             posend[index].append(self._pos)
             self._mo_types.append(index)
           except IndexError:
             pass
           try:
             index = 'Natural'
             self.find_string("MCSCF NATURAL ORBITALS")
             self.find_next_string("1          2")
             posend[index] = []
             posend[index].append(self._pos)
             self.find_next_string("LZ VALUE ANALYSIS FOR MOLECULAR ORBITALS")
             posend[index].append(self._pos)
             self._mo_types.append(index)
           except IndexError:
             pass
         self._mo_sets = {}
         for index  in self._mo_types:
            pos,end = posend[index]
            pos -= 1
            curatom = 1
            vectors = []
            while pos < end:
               pos += 1
               orbnum = self.text[pos].split()
               pos += 1
               line = self.text[pos].split()
               if len(line) == 0:
                  pos += 1
                  line = self.text[pos].split()
               try:
                 eigval = [ float(k) for k in line ]
                 pos += 1
                 syms = self.text[pos].split()
                 pos += 1
               except:
                 eigval = []
                 syms = []
                 pass
               line = self.text[pos].split()
               begin = pos
             # Build vectors
               pos = begin
               lmax = len(orbnum)+4
               for l in range(4,lmax):
                 pos = begin
                 line = self.text[pos]
                 line = line.replace('-',' -').split()
                 v = orbital()
                 v.set = index
                 v.basis = self.basis
                 if len(eigval)>0: v.eigenvalue = eigval[l-4]
                 if len(syms) >0: v.sym = syms[l-4]
                 while len(line) > 0 and pos < end:
                    try:
                       bid = int(line[2])
                    except ValueError:
                       line.insert(2,str(curatom))
                       self.text[pos] = ' '.join(line)
                    v.vector.append(float(line[l]))
                    pos += 1
                    line = self.text[pos]
                    line = line.replace('-',' -').split()
                 vectors.append(v)
            self._mo_sets[index] = vectors
         if 'BOYS' in self._mo_sets:
           self._mo_sets['BOYS'] += \
                self._mo_sets['RHF'][len(self._mo_sets['BOYS']):]

      return self._mo_sets

   def get_mulliken_mo(self):
      if self._mulliken_mo is None:
         posend = []
         vectors = []
         if self.spin_restrict:
           try:
              self.find_string("MULLIKEN ATOMIC POPULATION IN EACH MOLECULAR ORBITAL")
           except IndexError:
              return None
           pos = self._pos
           try:
              self.find_string("ATOMIC SPIN POPULATION")
           except IndexError:
              try:
                 self.find_next_string("POPULATIONS IN EACH AO")
              except IndexError:
                 return None
           end = self._pos
           pos += 2
           posend.append([pos,end])
         else:
           try:
              self.find_string("MULLIKEN ATOMIC POPULATION IN EACH MOLECULAR ORBITAL")
              self.find_next_string("ALPHA ORBITALS")
              pos = self._pos
              self.find_next_string("MULLIKEN ATOMIC POPULATION")
              end = self._pos
              pos += 2
              posend.append([pos,end])
              self.find_next_string("MULLIKEN ATOMIC POPULATION")
              self.find_next_string("BETA ORBITALS")
              pos = self._pos
              self.find_next_string("ATOMIC SPIN POP")
              end = self._pos
              pos += 2
              posend.append([pos,end])
           except IndexError:
              return None
         for pos, end in posend:
             while pos < end:
                try:
                  pos += 2
                  line = self.text[pos].split()
                  eigval = [ float(k) for k in line ]
                  pos += 2
                  line = self.text[pos].split()
                  begin = pos
                  lmax = len(eigval)+1
                  for l in range(1,lmax):
                    pos = begin
                    line = self.text[pos].split()
                    v = orbital()
                    v.set = "Mulliken_MO"
                    v.basis = self.basis
                    v.eigenvalue = eigval[l-1]
                    while len(line) > 1 and pos < end:
                       v.vector.append(float(line[l]))
                       pos += 1
                       line = self.text[pos].split()
                    vectors.append(v)
                    pos += 1
                except:
                  return None
         self._mulliken_mo = vectors
      return self._mulliken_mo

   def get_mulliken_ao(self):
       if self._mulliken_ao is None:
          try:
             self.find_string(" POPULATIONS IN EACH AO")
          except IndexError:
             return None
          self._pos += 2
          pos = self._pos
          self.find_next_string("----")
          end = self._pos - 1
          line = self.text[pos][26:].split()
          vm = orbital()
          vl = orbital()
          vm.set = "Mulliken_AO"
          vl.set = "Lowdin_AO"
#         vm.basis = self.basis
#         vl.basis = self.basis
          vm.eigenvalue = 0.
          vl.eigenvalue = 0.
          while pos < end:
             value = float(line[0])
             vm.vector.append(value)
             vm.eigenvalue += value
             value = float(line[1])
             vl.vector.append(value)
             vl.eigenvalue += value
             pos += 1
             line = self.text[pos][26:].split()
          self._mulliken_ao = vm
          self._lowdin_ao = vl
       return self._mulliken_ao

   def get_lowdin_ao(self):
       if self._lowdin_ao is None:
          self.get_mulliken_ao()
       return self._lowdin_ao


   def get_mulliken_atom(self):
       if self._mulliken_atom is None:
          try:
             self.find_string(" TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS")
          except IndexError:
             return None
          self._pos += 2
          pos = self._pos
          self.find_next_string("-----------")
          end = self._pos - 1
          line = self.text[pos].split()
          vm = orbital()
          vl = orbital()
          vm.set = "Mulliken_Atom"
          vl.set = "Lowdin_Atom"
#         vm.basis = self.geometry
#         vl.basis = self.geometry
          vm.eigenvalue = 0.
          vl.eigenvalue = 0.
          while pos < end:
             value = float(line[2])
             vm.vector.append(value)
             vm.eigenvalue += value
             value = float(line[4])
             vl.vector.append(value)
             vl.eigenvalue += value
             pos += 1
             line = self.text[pos].split()
          self._mulliken_atom = vm
          self._lowdin_atom = vl
       return self._mulliken_atom

   def get_lowdin_atom(self):
       if self._lowdin_atom is None:
          self.get_mulliken_atom()
       return self._lowdin_atom

   def get_two_e_int_ao(self):
      if self._two_e_int_ao is None:
         try:
            self.find_string("2 ELECTRON INTEGRALS")
         except IndexError:
            return None
         self.find_next_string("BYTES/INTEGRAL")
         pos = self._pos
         line = self.text[pos].split()
         Nint = int( line[1] )
         Nbytes = int( line[6] )
         assert ( Nbytes == 12 or Nbytes == 16 )
         try:
           file = gamessFile.two_e_int_array(self.two_e_int_filename)
         except IOError:
           return None
         if Nbytes == 12:
            tc = 1
            tci = 'b'
         elif Nbytes == 16:
            tc = 2
            tci = 'h'
         file.form = 'l'+4*Nint*tci+Nint*'d'
         self._two_e_int_ao = file
      #TODO
      return self._two_e_int_ao

   def get_num_alpha(self):
      if self._num_alpha is None:
         self._num_alpha = (self.num_elec + self.multiplicity-1)/2
      return self._num_alpha 

   def get_num_beta(self):
      if self._num_beta is None:
         self._num_beta = (self.num_elec - self.multiplicity+1)/2
      return self._num_beta 

   def get_determinants_mo_type(self):
      if self._determinants_mo_type is None:
        # Mono-determinant case
        if self.options['SCFTYP'] == 'RHF' or\
           self.options['SCFTYP'] == 'UHF' or\
           self.options['SCFTYP'] == 'ROHF':
           self._determinants_mo_type = self.mo_types[-1]
        # Multi-determinant case TODO
        elif self.options['SCFTYP'] == 'MCSCF':
           self._determinants_mo_type = 'MCSCF'
        # Multi-determinant case TODO
        elif self.options['CITYP'] != 'NONE':
           self._determinants_mo_type = self.mo_types[-1]
      return self._determinants_mo_type

   def get_csf_mo_type(self):
      if self._csf_mo_type is None:
         self._csf_mo_type = self.determinant_mo_type
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
          self._determinants = determinants
      return self._determinants

   def get_csf(self):
      if self._csf is None:

        # Mono-determinant case
        csf = []
        if self.options['SCFTYP'] == 'RHF' or\
           self.options['SCFTYP'] == 'UHF' or\
           self.options['SCFTYP'] == 'ROHF':
           new_csf = CSF()
           new_spin_a = []
           new_spin_b = []
           for i in range(self.num_alpha):
               new_spin_a.append(i)
           for i in range(self.num_beta):
               new_spin_b.append(i)
           new_csf.append(1.,new_spin_a,new_spin_b)
           csf.append(new_csf)
           self._determinants_mo_type = self.mo_types[-1]

        # Multi-determinant case
        if self.options['CITYP'] == 'GUGA' or self.options['SCFTYP'] == 'MCSCF':
           try:
             self.find_string('SYMMETRIES FOR THE')
             pos = self._pos
             buffer = self.text[pos].split()
             ncore = int (buffer[3])
             nact  = int (buffer[5])
             self.find_string('TOTAL NUMBER OF INTEGRALS')
             end = self._pos
             self.find_string('DETERMINANT CONTRIBUTION TO CSF')
             self.find_next_string('CASE VECTOR')
             pos = self.pos
             while pos < end:
               pos += 4
               if len (self.text[pos].split()) > 0:
                 this_csf = CSF()
                 while not self.text[pos].startswith(' CASE VECTOR') and len(self.text[pos])>1:
                   buffer = self.text[pos].split('=')[1].split(':')
                   mostr = buffer[1].replace("-"," -").split()
                   mo = []
                   for i in mostr:
                     mo.append (int(i))
                   p_count=0
                   for j in range ( len(mo) ):
                     if mo[j] < 0:
                       for l in range ( j+1, len(mo) ):
                         if mo[l] > 0:
                           p_count+=1
                   coef = float (buffer[0])*(-1.0)**p_count
                   tempcsf_a = []
                   tempcsf_b = []
                   for i in range(ncore):
                     tempcsf_a.append(i)
                     tempcsf_b.append(i)
                   for i in range(len(mo)):
                     orb = int (mo [i])
                     if orb > 0:
                       tempcsf_a.append(orb-1)
                     else:
                       tempcsf_b.append(-orb-1)
                   this_csf.append(coef,tempcsf_a,tempcsf_b)
                   pos+=1
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
        if self.options['SCFTYP'] == 'RHF' or\
           self.options['SCFTYP'] == 'UHF' or\
           self.options['SCFTYP'] == 'ROHF':
           self._det_coefficients = [ [1.] ]
        # Multi-determinant case
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
        if self.options['SCFTYP'] == 'RHF' or\
           self.options['SCFTYP'] == 'UHF' or\
           self.options['SCFTYP'] == 'ROHF':
           self._csf_coefficients = [ [1.] ]
        # Multi-determinant case 
        elif self.options['CITYP'] == 'GUGA' or self.options['SCFTYP'] == 'MCSCF':
           self._csf_coefficients = [ [] for i in range(self.num_states) ]
           self._pos = 0
           for state in range(self.num_states):
             cur_csf=0
             self.find_next_string("STATE #")
             self._pos += 4
             pos = self._pos
             try:
               while True:
                 line = self.text[pos].split()
                 csf=int(line[0])
                 for c in range(cur_csf,csf-1):
                    self._csf_coefficients[state].append(0.)
                 cur_csf=csf
                 self._csf_coefficients[state].append(float(line[1]))
                 pos += 1
             except:
               pass

      return self._csf_coefficients

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


   def get_pseudo(self):
     if self._pseudo is None:
        try:
          self.find_string("ECP POTENTIALS")
        except IndexError:
          return None
        pos = self._pos
        try:
          self.find_string("THE ECP RUN REMOVES")
        except IndexError:
          return None
        end = self._pos
        pos += 3
        pseudo_read = []
        while pos < end:
          line = self.text[pos][35:].split()
          pos += 1
          ecp = {}
          try:
            atom = line[0]
          except:
            continue
          try:
            ecp["zcore"] = int(line[3])
            ecp["atom"]  = atom
          except ValueError:  # Same as ...
            print line
            ecp = dict( pseudo_read[ int(line[6])-1 ] )
            ecp["atom"]  = atom
          else:
            lmax  = int(line[6])
            ecp["lmax"] = lmax
            for l in range(lmax+1):
              line = self.text[pos].split()
              l = int(line[2])
              pos += 1
              contraction = []
              while True:
                line = self.text[pos].split()
                try:
                  i = int(line[0])
                except ValueError:
                  break
                except IndexError:
                  break
                coef = float(line[1])
                n    = int(line[2])
                zeta = float(line[3])
                contraction.append( (coef, n, zeta) )
                pos += 1
              ecp[str(l)] = contraction
          pseudo_read.append(ecp)

        self._pseudo = pseudo_read
     return self._pseudo



# Properties
# ----------
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

   atom_to_ao_range = property(fget=get_atom_to_ao_range)

   for i in "two_e_int_mo_filename two_e_int_ao_filename".split():
     exec """
def get_%(i)s(self): return self._%(i)s
def set_%(i)s(self,value): self._%(i)s = value
%(i)s = property(fget=get_%(i)s,fset=set_%(i)s) """%locals()


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

def gamess_write_contrl(res,file, runtyp="ENERGY",guess="HUCKEL",
   exetyp="RUN",scftyp=None,moread=False,NoSym=False,units="BOHR"):
   print >>file, " $CONTRL"
   print >>file, "   EXETYP=",exetyp
   print >>file, "   COORD= UNIQUE",  " UNITS=",units
   print >>file, "   RUNTYP=", runtyp
   if NoSym:
     print >>file, "   NOSYM=1"
   if scftyp is None:
     if res.num_alpha == res.num_beta:
       scftyp="RHF"
     else:
       scftyp="ROHF"
   print >>file, "   SCFTYP=", scftyp
   print >>file, "   MULT=", res.multiplicity
   print >>file, "   ICHARG=",str(int(res.charge))
   print >>file, " $END"
   print >>file, ""
   print >>file, " $GUESS"
   if moread:
     print >>file, "  GUESS=MOREAD"
     print >>file, "  NORB=",len(res.mo_sets[res.mo_types[-1]])
   else:
     print >>file, "  GUESS=HUCKEL"
   print >>file, " $END"


def gamess_write_vec(res, file):
#  if motype is None:
#    motypelist = res.mo_types
#  else:
#    motypelist = [motype]
#  for motype in motypelist:
     motype = res.mo_types[-1]
     print >>file, motype, "MOs"
     print >>file, " $VEC"
     moset = res.mo_sets[motype]
     moindex = res.closed_mos + res.active_mos + \
               res.virtual_mos
     for idx,moi in enumerate(moindex):
       mo = moset[moi]
       v = list(mo.vector)
       line = 1
       while (len(v) > 0):
         fields = []
         monum = "%4d"%(idx+1 )
         monum = monum[-2:]
         fields.append ( monum )
         linum = "%4d"%(line)
         linum = linum[-3:] 
         fields.append ( linum )
         for n in range(5):
           try:
             fields.append ( "%15.8E"%v.pop(0) )
           except IndexError:
             pass
         print >>file, ''.join(fields)
         line += 1
     print >>file, " $END"

def gamess_write_data(res, file):
  print >>file, " $DATA"
  print >>file, res.title
  try:
    point_group = res.point_group
    N = point_group[1]
    if N != '1':
     point_group = point_group[0]+'N'+point_group[2:].capitalize()
     point_group += '    '+N+'\n'
  except TypeError:
    point_group = 'C1'
    N = '1'
  print >>file, point_group
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
      elif sym == 'd+0' :
         sym = 'D'
         doPrint = True
      elif sym == 'xx' :
         sym = 'D'
         doPrint = True
      elif sym == 'f+0' :
         sym = 'F'
         doPrint = True
      elif sym == 'xxx' :
         sym = 'F'
         doPrint = True
      elif sym == 'g+0' :
         sym = 'G'
         doPrint = True
      elif sym == 'xxxx' :
         sym = 'G'
         doPrint = True
      elif sym == 'xxxxx' :
         sym = 'H'
         doPrint = True
      elif sym == 'xxxxxx' :
         sym = 'I'
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

def gamess_write_integrals(res,file):
  filename = file.name.replace('.inp','.moints')
      
def write_gamess(file,stdout):
  gamess_write_contrl(file,stdout)
  x = file.basis
  x = file.geometry
  gamess_write_data(file,stdout)
  gamess_write_vec(file,stdout)

  
fileTypes.append(gamessFile)

if __name__ == '__main__':
   main(gamessFile)
