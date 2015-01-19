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

molproFile_defined_vars = [ "date", "version", "machine", "memory", "disk",\
                "cpu_time", "author", "title", "units", "methods", \
                "spin_restrict", "conv_threshs", "energies", \
                "one_e_energies", "two_e_energies", \
                "kin_energies", "virials", "point_group", "num_elec", \
                "charge", "multiplicity","nuclear_energy","dipole","geometry",\
                "basis","mo_sets","mo_types",\
                "determinants", "num_alpha", "num_beta",\
                "closed_mos", "active_mos", "virtual_mos", \
                "determinants_mo_type", "det_coefficients", \
                "csf_mo_type", "csf_coefficients", "symmetries", "occ_num", \
                "csf", "num_states", "num_orb_sym" ]

methods_with_orbitals = ['RHF-SCF', 'MULTI']

class molproFile(resultsFile):
   """ Class defining the molpro file.
   """

   local_vars = list(local_vars)
   defined_vars = list(molproFile_defined_vars)

   exec get_data('date',"DATE: ",'-3:-2') in locals()
   exec get_data('point_group',"Point group",'2:') in locals()
   exec get_data('version',"    Version",'1:2') in locals()
   exec get_data('nuclear_energy',"NUCLEAR REPULSION ENERGY", '3:4',"float") in locals()

   def get_num_elec(self):
     if self._num_elec is None:
       self._num_elec = self.num_alpha + self.num_beta
     return self._num_elec 

   def get_charge(self):
     if self._charge is None:
       self._charge = 0.
       for a in self.geometry:
         self._charge += a.charge
       self._charge -= self.num_elec
     return self._charge

   def get_multiplicity(self):
     if self._multiplicity is None:
       self._multiplicity = self.num_alpha - self.num_beta + 1
     return self._multiplicity

   def get_author(self):
     if self._author is None:
       buffer = ' '.join(self.text[-10:])
       regexp = re.compile(r"user=(.+),")
       matches = regexp.findall(buffer)
       if matches != []:
         self._author = matches[0]
     return self._author

   def get_machine(self):
     if self._machine is None:
       regexp = re.compile(r"[^/]+/([^\(]+).+DATE:")
       for line in self.text:
         matches = regexp.findall(line)
         if matches != []:
           self._machine = matches[0]
           return self._machine
     return self._machine

   def get_disk(self):
     if self._disk is None:
       buffer = ' '.join(self.text[-100:])
       regexp = re.compile(r"DISK USED +\* +([0-9]+\.[0-9]* +[^ \n]+)")
       matches = regexp.findall(buffer)
       try:
         self._disk = matches[-1].title()
       except IndexError:
         return None
     return self._disk

   def get_cpu_time(self):
     if self._cpu_time is None:
       buffer = ' '.join(self.text[-100:])
       regexp = re.compile(r"CPU TIMES +\* +([0-9]+\.[0-9]*) +[^ \n]+")
       matches = regexp.findall(buffer)
       try:
         self._cpu_time = matches[-1] + " s"
       except IndexError:
         return None
     return self._cpu_time

   def get_memory(self):
     if self._memory is None:
       buffer = ' '.join(self.text[-100:])
       regexp = re.compile(r"SF USED +\* +([0-9]+\.[0-9]* +[^ \n]+)")
       matches = regexp.findall(buffer)
       try:
         self._memory = matches[-1].title()
       except IndexError:
         return None
     return self._memory

   def get_title(self):
     if self._title is None:
       try:
         self.find_string("LABEL")
       except IndexError:
         return None
       pos = self._pos
       self._title = ' '.join(self.text[pos].split()[2:]).strip()
     return self._title

   def get_symmetries(self):
     if self._symmetries is None:
       try:
         self.find_string("NUMBER OF CONTRACTIONS")
       except IndexError:
         return None
       pos = self._pos
       regexp = re.compile(r" +([0-9]+)([A-Z][^ ]*) +")
       self._symmetries = [ [j,int(i)] for i,j in regexp.findall(self.text[pos]) ]
     return self._symmetries

   def get_units(self):
     if self._units is None:
       self._units = "BOHR"
     return self._units

   def get_methods(self):
     if self._methods is None:
        regexp = re.compile(r"^1PROGRAM *\* *([^ ]+)")
        methods = []
        for line in self.text:
          program = regexp.findall(line)
          if program != []:
           if program[0] in methods_with_orbitals:
            methods += program
        self._methods = methods
     return self._methods

   def get_occ_num(self):
     if self._occ_num is None:
        occ = {}
        program = []
        progFound=False
        doRead=False
        regexp  = re.compile(r"^1PROGRAM *\* *([^ ]+)")
        regexp2 = re.compile(r"^ +Orb +Occ +")
        regexp3 = re.compile(r"^ +[0-9]*\.[0-9]+ +([0-2]*\.[0-9]+|[0-2]|\+|\-) +[^ ]")
        for line in self.text:
          if not progFound:
            program = regexp.findall(line)
            if program != []:
             if program[0] in methods_with_orbitals:
              while program[0] in occ.keys():
                program[0] += 'x'
              occ[program[0]] = []
              progFound=True
          else:
            buffer = line[:25]
            if doRead:
              if buffer.strip() != "":
                matches = regexp3.findall(buffer)
                if matches != []:
                  if matches[0] == '+':
                    res = 1.
                  else:
                    res = float(matches[0])
                  occ[program[0]].append(res)
                else:
                  doRead = False 
                  progFound = False 
            elif regexp2.match(buffer):
              doRead = True
        if occ != {}:
          self._occ_num = occ
     return self._occ_num 

   def get_spin_restrict(self):
       if self._spin_restrict is None:
          self._spin_restrict = True
          for m in self.methods:
             if "U" in m:
                self._spin_restrict = False
       return self._spin_restrict


   def get_conv_threshs(self):
       if self._conv_threshs is None:
          self._conv_threshs = []
          regexp = \
            re.compile(r" conv.* thre.* +.*([0-9]+\.[0-9]*[^ ]*) *\(energy\)",
                            re.I)
          for line in self.text:
            matches = regexp.findall(line)
            if matches != []:
               self._conv_threshs += [float(matches[0])]
       return self._conv_threshs

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
       pos = len(self.text)-1
       while not self.text[pos].startswith(" PROGRAMS") and pos>0:
          pos -=1
       while not self.text[pos].startswith(" ****"):
          pos += 1
       pos +=1 
       regexp = re.compile(r"^\s+([|\-][0-9]+\.[0-9]+\s*)+")
       found = False
       while not found:
          result = regexp.match(self.text[pos])
          if result != None: 
             found = True
             buffer = self.text[pos].split()
             for e in buffer:
               self._energies.insert(0,float(e))
          pos += 1
       program = []
       progFound=False
       regexp  = re.compile(r"^1PROGRAM *\* *([^ ]+)")
       regexp2 = re.compile(r"^ (Virial|One.elec|Two.elec|Kinetic)[^ ]* +[^ ]+ +(.+)")
       progindex = -1
       for line in self.text:
         if not progFound:
           program = regexp.findall(line)
           if program != []:
            if program[0] in methods_with_orbitals:
             progindex += 1
             progFound=True
         else:
           try:
             quantity, value = regexp2.findall(line)[0]
             if quantity == "Virial":
                quantity = self._virials
                progFound=False
             elif quantity[0:3] == "One":
                quantity = self._one_e_energies
             elif quantity[0:3] == "Two":
                quantity = self._two_e_energies
             elif quantity == "Kinetic":
                quantity = self._kin_energies
             while len(quantity) <= progindex:
                quantity.append(None)
             quantity[progindex] = float(value)
           except:
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

   def get_kin_energies(self):
     if self._kin_energies is None:
        self.get_energies()
     return self._kin_energies

   def get_virials(self):
     if self._virials is None:
        self.get_energies()
     return self._virials

   def get_num_alpha(self):
     if self._num_alpha is None:
        try:
           self.find_string(" NUMBER OF ELECTRONS")
           pos = self._pos
           buffer = self.text[pos].split()
           self._num_alpha = int(buffer[3][:-1])
           self._num_beta = int(buffer[4][:-1])
        except IndexError:
           pass
     return self._num_alpha

   def get_num_beta(self):
     if self._num_beta is None:
       self.get_num_alpha()
     return self._num_beta

   def get_geometry(self):
     if self._geometry is None:
        try:
           self.find_string(" ATOMIC COORDINATES")
           pos = self._pos
           pos += 4
           self._geometry = []
           buffer = self.text[pos].split()
           while len(buffer) > 1:
              temp = atom()
              temp.name   = buffer[1]
              temp.charge = float(buffer[2])
              temp.coord  = (float(buffer[3]), float(buffer[4]), float(buffer[5]))
              temp.basis  = []
              self._geometry.append(temp)
              pos += 1
              buffer = self.text[pos].split()
        except IndexError:
           pass
     return self._geometry

   def get_dipole(self):
     if self._dipole is None:
        try:
           self.find_last_string("DIPOLE MOMENTS:")
           pos = self._pos
           self._dipole = []
           buffer = self.text[pos].split(':')[1]
           for d in buffer.split()[:3]:
             self._dipole.append( float(d) )
        except IndexError:
           pass
     return self._dipole
   
   def get_num_orb_sym(self):
     try:
       if self._num_orb_sym is None:
         self.get_basis()
     except:
       self._num_orb_sym = None
       self.get_basis()
     return self._num_orb_sym

   def get_basis(self):
     if self._basis is None:
        try:
           self.find_string("BASIS DATA")
        except IndexError:
           return None
        pos = self._pos
        try:
           self.find_next_string("NUCLEAR CHARGE")
        except IndexError:
           return None
        end = self._pos-1
        num_orb_sym = []
        for k in self.symmetries:
           num_orb_sym.append([0,0])
        pos += 4
        basis_read = []
        index, sym, nuc, type, expo, coef = "", "", "", "", "", ""
        basis = []
        currfunc = 0
        salcao2mo = []
        contrmax=0
        buffer = self.text[pos:end+1]
        text = []
        for line in buffer:
          if line[:41].strip() == "":
            text[-1] = text[-1][:-1] + line[41:]
          else:
            text.append(line)
        pos=0
        end = len(text)
        text.append(80*" ")
        while pos < end:
           begin = pos
           line = text[begin]
           nfunc = len( line[41:].split() )
           nuc = []
           try:
              idx4 = line.find('.')
              idx5 = idx4+1
              idx6 = idx4+2
              currfunc = int(line[:idx4])
              geom = int( line[idx5:idx6] )
           except:
             line = str(currfunc)+'.'+str(geom)+line[8:]
           end2 = pos
           while end2<end and int( line[:idx4] ) < currfunc+nfunc and \
             int( line[idx5:idx6] ) == geom:
             if  line[11:15].strip() != '':
               nuc.append( (int(line[11:15])-1, line[15:22].strip() ) )
             end2 += 1
             line = text[end2]
             if line[:idx4] == '    ':
               line = "%4d.%d  "%(currfunc,geom)+line[8:]
           line = text[begin]

           
           contr = []
           for k in range(nfunc):
             for NucType in nuc:
               c = contraction()
               contr.append(c)

           thiscontr=0
           if nfunc > contrmax: contrmax = nfunc
           for k in range(nfunc):
             smo = []
             num_orb_sym[geom-1][0] += 1
             num_orb_sym[geom-1][1] += len(nuc)
             for NucType in nuc:
               c = contr[thiscontr]
               thiscontr += 1
               sign=1.
               pos = begin
               line = text[pos]
               atom = self.geometry[NucType[0]]
               type = str(NucType[1])
               if type[0] == '-':
                 type = type[1:]
                 sign = -1.
               if type[0].isdigit():
                  type = type[1:]
               try:
                 if not type[1].isdigit() :
                   if type[0].lower() in "pdfghi": type = type[1:]
               except:
                 if type[0].lower() in "pdfghi": type = type[1:]
               smo.append(sign/sqrt(float(len(nuc))))
               c.center = atom.coord
               type = normalize_basis_name(type)
               c.sym = type
               while pos<end2:
                 expo = line[24:40].strip()
                 coef = line[41:].split()
                 if expo == '':
                   pos = end2
                   break
                 gauss = gaussian()
                 gauss.center = atom.coord
                 gauss.expo = float(expo)
                 gauss.sym  = type
                 if float(coef[k]) != 0.:
                   c.append(float(coef[k]),gauss)
                 pos += 1
                 line = text[pos]
             salcao2mo.append( smo )
           for c in contr:
  #          c.normalize()
             basis.append(c)
        self._num_orb_sym = num_orb_sym
        self._basis = basis
        self.salcao2mo = salcao2mo
        for f in basis:
          for at in self.geometry:
            if f.center is at.coord:
              at.basis.append(f)
        if (len(self.symmetries)>1) and contrmax > 1:
          file = open('basis.molpro','w')
          molpro_write_input(self,file)
          file.close()
          print >>sys.stderr, """
Warning: there is a bug in the molpro output. Run the 'basis.molpro'
input file with molpro to generate a correct output for the basis set.
"""
     return self._basis

   def get_mo_types(self):
      if self._mo_types is None:
        self.get_mo_sets()
      return self._mo_types 


   def get_mo_sets(self):
      if self._mo_sets is None:
        global methods_with_orbitals
        self._mo_sets = {}
        self._mo_types = []
        doRead=False
        regexp  = re.compile(r"^1PROGRAM *\* *([^ ]+)")
        regexp2 = re.compile(r"^ +Orb +Occ +")
        rstring  = r"^ {,4}([0-9]+\.[0-9]+) +"             # Orb
        rstring += r"([0-2| ]\.[0-9]{,5}|[0-2]|\+|\-) *" # Occ
        rstring += r"(-*[0-9|-]*\.[0-9]+) +" # Energy
        regexp3 = re.compile(rstring)
        rstring = r"(([0-9|\-]*\.[0-9]{6} *)+)" # Coefs
        regexp4 = re.compile(rstring)
        rstring = r"(^ {30,}([0-9] [123456789xyzspdfghijk\+\-]+ *)+)" # Labels
        regexp5 = re.compile(rstring)
        rstring = r"(^ \*{10,}|orbital dump)" # *********
        regexp6 = re.compile(rstring)
        rstring = " +[0-9|-]*\.[0-9]+\.[0-9]*"
        regexp7 = re.compile(rstring)
        text = list(self.text)
        curr_sym = self.basis[0].sym
        v=None
        ReadingLabels = False
        PreviousBlank = False
        for iline, line in enumerate(text):
          buffer = regexp7.findall(line)
          while buffer != []:  # Insert spaces in "7851.548715699.7413"
             if ":" in line:
               break
             begin = line.find(buffer[0])
             end = begin+len(buffer[0])
             line = line[0:end-10]+" "+line[end-10:]
             buffer = regexp7.findall(line)
          line = line[:11]+line[11:32].replace ("**********"," 999999.9999") \
                          +line[32:].replace('-',' -')
          # Find method
          program = regexp.findall(line)
          if program != []:
           if program[0] in methods_with_orbitals:
            add_method = False
            while program[0] in self._mo_sets:
              program[0] += 'x'
              add_method = True
            if add_method and program[0] not in methods_with_orbitals:
              methods_with_orbitals += [program[0]]
            index = program[0]
            self._mo_types.append(index)
            self._mo_sets[index] = []
          if doRead:
            # Find labels
            buffer = regexp6.findall(line)
            if buffer != []:
              doRead=False
            elif line.strip() == "":
              if symcount == len(self.symmetries):
                if PreviousBlank:
                  doRead = False
                  PreviousBlank = False
                else:
                  if not ReadingLabels:
                    PreviousBlank = True
              if v is not None:
                for l,s in enumerate(self.num_orb_sym):
                   if l > symcount-1 :
                     for k in range(s[1]):
                       v.vector.append(0.)
                self._mo_sets[index].append(v)
                v = None
            elif line.strip() != "":
              buffer = regexp5.findall(line)
              PreviousBlank = False
              if buffer != []:
                if not ReadingLabels:
                  ReadingLabels = True
              else:
                ReadingLabels = False
              # Find "Orb Occ Energy" values to initiate a new orbital
              buffer = regexp3.findall(line)
              if buffer != []:
                buffer = buffer[0]
                symcount = int(buffer[0].split('.')[1])
                bf = 0
                v = orbital()
                v.set = index
                v.basis = self.basis
                v.sym = self.symmetries[ int(buffer[0].split('.')[1])-1 ][0]
                v.eigenvalue = float(buffer[2])
                if v.eigenvalue == 999999.9999:
                  print "Warning line", iline+1, ": '**********' in orbital eigenvalues"
                for l,s in enumerate(self.num_orb_sym):
                  if l < symcount-1 :
                    bf += s[0]
                    for k in range(0,s[1]):
                      v.vector.append(0.)
              # Find coefficient values to continue the orbital
              buffer = regexp7.findall(line[30:])
              while buffer != []:  # Insert spaces in "7851.548715699.7413"
                 begin = line[0:].find(buffer[0])
                 end = begin+len(buffer[0])
                 line = line[0:end-10]+" "+line[end-10:]
                 buffer = regexp7.findall(line)
              buffer = regexp4.findall(line[30:])
              if buffer != []:
                for x in buffer[0][0].split():
                  for k in self.salcao2mo[bf]:
                    try:
                      x2 = k*float(x)
                    except ValueError:
                      print "Error line", iline+1, ": orbital coefficients"
                      sys.exit(1)
                    v.vector.append(x2)
                  bf+=1
          elif regexp2.match(line):
            if self._mo_sets[index] == []:
              doRead = True
            symcount=0
      return self._mo_sets

   def get_determinants_mo_type(self):
      if self._determinants_mo_type is None:
        self._determinants_mo_type = self.mo_types[-1]
      return self._determinants_mo_type

   def get_csf_mo_type(self):
      if self._csf_mo_type is None:
        self._csf_mo_type = self.determinants_mo_type
      return self._determinants_mo_type

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
        csfcoef = [[]]
        index = len(self.methods)-1
        while self.mo_types[index] not in methods_with_orbitals and index >0:
           index -= 1
        # Mono-determinant case
        if self.mo_types[index].startswith("RHF-SCF"):
           new_csf = CSF()
           new_spin_det_a = []
           new_spin_det_b = []
           occ = self.occ_num[self.determinants_mo_type]
           for i,v in enumerate(occ):
                   if v == 1.:
                           new_spin_det_a.append(i)
                   elif v == 2.:
                           new_spin_det_a.append(i)
                           new_spin_det_b.append(i)
                   else:
                           assert v == 0.
           new_csf.append(1.,new_spin_det_a,new_spin_det_b)
           csf.append(new_csf)
           csfcoef[0].append(1.)
           self._determinants_mo_type = self.mo_types[index]

        # Multi-determinant case
        elif self.mo_types[index].startswith("MULTI"):
             try:
               self.find_last_string('Number of closed-shell orbitals')
               pos = self._pos
               buffer = self.text[pos].split()
               ncore0 = [ int(a) for a in buffer[6:-1] ]
             except IndexError:
               ncore0 = [ 0 for a in self.symmetries ]

             self.find_last_string('Number of active')
             pos = self._pos
             buffer = self.text[pos].split()
             nact0  = [ int(a) for a in buffer[6:-1] ]

             self.find_next_string('Number of external')
             mo_count = [ 0 for i in self.num_orb_sym ] 
             j = 0
#            curr_sym = self.mo_sets[self.determinants_mo_type][0].sym[0]
             curr_sym = self.mo_sets[self.determinants_mo_type][0].sym
             for i in self.mo_sets[self.determinants_mo_type]:
              if i.sym is not None:
               if i.sym != curr_sym:
                 curr_sym = i.sym
#              if i.sym[0] != curr_sym:
#                curr_sym = i.sym[0]
                 j+=1
               mo_count[j] += 1
             next0  = [ mo_count[i] - nact0[i] - ncore0[i] for i in range(len(self.num_orb_sym)) ]
             self.find_next_string('TOTAL ENERGIES')
             end = self._pos -2
             self.find_prev_string('CI vector')
             pos = self._pos +2
             while pos < end:
              pos += 1
              this_csf = CSF()
              tempcsf_a = []
              tempcsf_b = []
              buffer = [ ''.join(self.text[pos].split()[:-1]),self.text[pos].split()[-1]]
              mo, coef = buffer
              coef = float(coef)
              ncoreold=0
              nactold=0
              mo2=""
              for ((ncore,nact),next) in zip(zip(ncore0,nact0),next0):
               for i in range(ncore):
                 tempcsf_a.append(ncoreold+i)
                 tempcsf_b.append(ncoreold+i)
                 mo2+="2"
               ncoreold += ncore
               for k,i in enumerate(mo[nactold:nactold+nact]):
                 if i in "2a+/":
                   tempcsf_a.append(ncoreold+k)
                 if i in "2b-\\":
                   tempcsf_b.append(ncoreold+k)
                 mo2 += i
               ncoreold += nact+next
               nactold  += nact

# Correction de la phase due au changement d'ordre des orbitales
              p_count=0
              for k,i in enumerate(mo2):
                if i in "2b-\\":
                  for j in mo2[k+1:]:
                    if j in "2a+/":
                      p_count+=1
              this_csf.append((-1.0)**p_count,tempcsf_a,tempcsf_b)

              csf.append(this_csf)
              csfcoef[0].append(coef)
        if csf != []:
           self._csf = csf
           self._csf_coefficients = csfcoef
      return self._csf


   def get_csf_coefficients(self):
      if self._csf_coefficients is None:
         self.get_csf()
      return self._csf_coefficients

   def get_closed_mos(self):
      if self._closed_mos is None:
         result = []
         occ = filter(lambda x: x>0., self.occ_num[self.determinants_mo_type])
         maxmo = len(occ)
         for orb in range(maxmo):
           present = True
           for det in self.determinants:
               if present:
                 for spin_det in det:
                   if orb not in det[spin_det]: present = False
               else:
                 break
           if present: result.append(orb)
         self._closed_mos = result
      return self._closed_mos

   def get_virtual_mos(self):
      if self._virtual_mos is None:
         result = []
         occ = filter(lambda x: x>0., self.occ_num[self.determinants_mo_type])
         maxmo = len(self.mo_sets[self.determinants_mo_type])
         for orb in range(maxmo):
           present = False
           for det in self.determinants:
             if not present:
               for spin_det in det:
                   if orb in det[spin_det]: present = True
             else:
               break
           if not present: result.append(orb)
           self._virtual_mos = result
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
        index = len(self.methods)-1
        while self.mo_types[index] not in methods_with_orbitals and index >0:
           index -= 1
        # Mono-determinant case
        if self.mo_types[index].startswith("RHF-SCF"):
           self._det_coefficients = [ [1.] ]
        # Multi-determinant case
        if self.mo_types[index].startswith("MULTI"):
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

   def get_num_states(self):
      if self._num_states is None:
        self._num_states=1
      return self._num_states

# Properties
# ----------
   exec build_property("num_orb_sym","Number of SALCAO/sym") in locals()
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

# Output Routines
# ---------------

def molpro_write_input(res,file):

  print >>file, " ***,", res.title
  geom = res.geometry
  print >>file, """
 print,basis;
 print,orbitals;
 gprint,civector;
 gthresh,printci=0.;

 geomtyp=xyz
 geometry={
 """, len(geom)
  print >>file, res.title
  i=1
  for at in geom:
    coord = []
    for x in at.coord: coord.append(x*a0)
    print >>file, "%10s  %15.8f %15.8f %15.8f" % tuple( \
     [(at.name+str(i)).ljust(10) ]+coord )
    i+=1
  print >>file, "}\n"
  print >>file, "basis={"
  lines = []
  for idx,at in enumerate(geom):
    # Find atom
    doPrint = True
    for contr in at.basis:
         line = ""
         # Find Sym
         sym = contr.sym
         if sym == 's' :
            sym = 's'
            doPrint = True
         elif sym == 'x' :
            sym = 'p'
            doPrint = True
         elif sym == 'd+0' :
            sym = 'd'
            doPrint = True
         elif sym == 'xx' :
            sym = 'd'
            doPrint = True
         elif sym == 'f+0' :
            sym = 'f'
            doPrint = True
         elif sym == 'xxx' :
            sym = 'f'
            doPrint = True
         elif sym == 'g+0' :
            sym = 'g'
            doPrint = True
         elif sym == 'xxxx' :
            sym = 'g'
            doPrint = True
         elif sym == 'h+0' :
            sym = 'h'
            doPrint = True
         elif sym == 'xxxxx' :
            sym = 'h'
            doPrint = True
         elif sym == 'i+0' :
            sym = 'i'
            doPrint = True
         elif sym == 'xxxxxx' :
            sym = 'i'
            doPrint = True
         elif sym == 'xxxxxxx' :
            sym = 'j'
            doPrint = True
         else:
            doPrint = False
         if doPrint:
           # Find exponents
           line += "%1s,%s%d"%(sym, at.name, idx+1)
           for i, p in enumerate(contr.prim):
              if not isinstance(p,gaussian):
                 raise TypeError("Basis function is not a gaussian")
              line += ",%f"%(float(p.expo))
           line += "\n"
           # Find coefs
           line += "c,1."+str(len(contr.prim))
           for c in contr.coef:
              line += ",%g"%(c)
           if line not in lines:
             lines.append(line)
  for line in lines:
    print >>file, line
  print >>file, "}\n",

  print >>file, "rhf;"
  print >>file, "---"
  
fileTypes.append(molproFile)

if __name__ == '__main__':
   main(molproFile)

