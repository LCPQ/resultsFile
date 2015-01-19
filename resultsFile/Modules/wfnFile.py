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

wfnFile_defined_vars = [ 
                "energies", "title", "units",\
                "virials", "num_elec", \
                "charge", "geometry",\
                "basis","mo_sets","mo_types", \
                "num_alpha", "num_beta",\
                "closed_mos", "active_mos", "virtual_mos", \
                "determinants",
                "determinants_mo_type", "det_coefficients", \
                "csf_mo_type", "csf_coefficients", "occ_num", \
                "csf", "num_states"] 

class wfnFile(resultsFile):
   """ Class defining the wfn file.
   """

   local_vars = list(local_vars)
   defined_vars = list(wfnFile_defined_vars)

   def get_title(self):
     if self._title is None:
       self._title = self.text[0].strip()
     return self._title

   def get_units(self):
     if self._units is None:
       self._units = "BOHR"
     return self._units

   def get_num_elec(self):
     if self._num_elec is None:
       self._num_elec = int(sum(self.occ_num["Wfn"]))
     return self._num_elec

   def get_occ_num(self):
     if self._occ_num is None:
        result = {}
        occ = []
        for line in self.text:
          if line.startswith("MO"):
            buffer = line.split()
            while buffer[0] != "OCC":
              buffer.pop(0)
            occ.append(float(buffer[3]))
        result["Wfn"] = occ
        self._occ_num = result
     return self._occ_num 

   def get_energies(self):
       if self._energies is None:
          buffer = self.text[-1]
          self._energies = [float(buffer.split("ENERGY")[1].split()[1])]
          self._virials =  [float(buffer.split("VIRIAL(-V/T)=")[1].split()[0])]
       return self._energies

   def get_virials(self):
     if self._virials is None:
        self.get_energies()
     return self._virials

   def get_charge(self):
     if self._charge is None:
        result = 0.
        for line in self.text:
          if "CHARGE" in line:
            result += float(line.split("CHARGE =")[1].strip())
        result -= self.num_elec
        self._charge = result
     return self._charge

   def get_geometry(self):
     if self._geometry is None:
        Natoms = int(self.text[1].split()[6])
        pos = 2
        begin = pos
        self._geometry = []
        buffer = self.text[pos].split()
        while pos<begin+Natoms:
           temp = atom()
           temp.name   = buffer[0]
           temp.charge = float(buffer[9])
           temp.coord  = (float(buffer[4]), float(buffer[5]), float(buffer[6]))
           temp.basis  = []
           self._geometry.append(temp)
           pos += 1
           buffer = self.text[pos].split()
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

   def get_basis(self):
     if self._basis is None:
        center = []
        type = []
        expo = []
        conversion = [ "S", "X", "Y", "Z" ]
        for line in self.text:
           if line.startswith('CENTRE'):
             buffer = line[20:].split()
             for b in buffer:
               center.append(int(b)-1)
           elif line.startswith('TYPE'):
             buffer = line[20:].split()
             for b in buffer:
               type.append(int(b))
           elif line.startswith('EXPONENTS'):
             buffer = line[10:].split()
             for b in buffer:
               expo.append(float(b.replace('D','e')))
        if 10 in type: # Cartesian d orbitals
          conversion += [ "XX", "YY", "ZZ", "XY", "XZ", "YZ" ]
        else:          # Spherical d
          conversion += [ "D2-", "D1-", "D0", "D1+", "D2+", "" ]
        if 20 in type: # Cartesian f orbitals
          conversion += [ "XXX", "YYY", "ZZZ", "YYX", "XXY", "XXZ",
                          "ZZX", "ZZY", "YYZ", "XYZ" ]
        else:          # Spherical f
          conversion += [ "F3-", "F2-", "F1-", "F0", "F1+", "F2+", 
                          "F3+", "", "", "" ]
        for i in range(len(type)):
          type[i] = conversion[type[i]-1]
        assert ( len(expo) == len(type) )
        assert ( len(expo) == len(center) )
        
        self._basis = []
        for i in range(len(expo)):
          contr = contraction()
          gauss = gaussian()
          atom = self.geometry[center[i]]
          gauss.center = atom.coord
          gauss.expo = expo[i]
          gauss.sym  = type[i]
          contr.append(1.,gauss)
          self._basis.append(contr)
     return self._basis

   def get_mo_types(self):
      if self._mo_types is None:
        self.get_mo_sets()
      return self._mo_types 

   def get_mo_sets(self):
      if self._mo_sets is None:
         self._mo_sets = {}
         index = 'Wfn'
         self._mo_types = [index]
         posend = {}
         k=0
         vectors = []
         for line in self.text:
           if line.startswith('MO'):
             if k>0:
               vectors.append(v)
             v = orbital()
             v.set = index
             v.basis = self.basis
             v.eigenvalue = float(line.split()[-1])
             v.sym = None
             k+=1
           elif line.startswith('END'):
             vectors.append(v)
             k=0
           elif k>0:
             for c in line.split():
               v.vector.append(float(c.replace('D','e')))
         basis = self.basis
         for v in vectors:
           for i in range(len(basis)):
             v.vector[i] *= basis[i].norm
         self._mo_sets[index] = vectors
      return self._mo_sets

   def get_num_alpha(self):
      if self._num_alpha is None:
         self._num_alpha = self.num_elec-self.num_beta
      return self._num_alpha 

   def get_num_beta(self):
      if self._num_beta is None:
         self._num_beta = self.num_elec/2
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
        new_csf = CSF()
        new_spin_det_a = []
        new_spin_det_b = []
        ntot = 0
        for i,n in enumerate(self.occ_num['Wfn']):
          if n == 1.0:
            if ntot < self.num_alpha:
              new_spin_det_a.append(i)
            else:
              new_spin_det_b.append(i)
            ntot += 1
          elif n == 2.0:
            new_spin_det_a.append(i)
            new_spin_det_b.append(i)
        new_csf.append(1.,new_spin_det_a,new_spin_det_b)
        csf.append(new_csf)
        self._determinants_mo_type = self.mo_types[-1]
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
        self._det_coefficients = [ [1.] ]
      return self._det_coefficients

   def get_csf_coefficients(self):
      if self._csf_coefficients is None:
        self._csf_coefficients = [ [1.] ]
      return self._csf_coefficients

   def get_num_states(self):
      if self._num_states is None:
        self._num_states=1
      return self._num_states

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


# Output Routines
# ---------------

import string

basis_correspond = {}
basis_correspond['s'  ] = 1
basis_correspond['x'  ] = 2
basis_correspond['y'  ] = 3
basis_correspond['z'  ] = 4
basis_correspond['xx' ] = 5
basis_correspond['yy' ] = 6
basis_correspond['zz' ] = 7
basis_correspond['xy' ] = 8
basis_correspond['xz' ] = 9
basis_correspond['yz' ] = 10
basis_correspond['xxx'] = 11
basis_correspond['yyy'] = 12
basis_correspond['zzz'] = 13
basis_correspond['xyy'] = 14
basis_correspond['xxy'] = 15
basis_correspond['xxz'] = 16
basis_correspond['xzz'] = 17
basis_correspond['yzz'] = 18
basis_correspond['yyz'] = 19
basis_correspond['xyz'] = 20
basis_correspond['d-2'] = 5
basis_correspond['d-1'] = 6
basis_correspond['d+0'] = 7
basis_correspond['d+1'] = 8
basis_correspond['d+2'] = 9
basis_correspond['f-3'] = 11
basis_correspond['f-2'] = 12
basis_correspond['f-1'] = 13
basis_correspond['f+0'] = 14
basis_correspond['f+1'] = 15
basis_correspond['f+2'] = 16
basis_correspond['f+3'] = 17

def wfn_write(res,file,MO_type=None):
  print >>file, " "+res.title.strip()
  print >>file, "GAUSSIAN %14d MOL ORBITALS %6d PRIMITIVES %8d NUCLEI" \
    %(len(res.closed_mos+res.active_mos), len(res.uncontracted_basis), len(res.geometry)),

  # write geom
  geom  = res.geometry
  f=1.
  if res.units == 'ANGS':
    f=1./a0
  for i,atom in enumerate(geom):
   name = atom.name
   for d in string.digits:
     name = name.replace(d,'')
   print >>file, "\n%3s %4d    (CENTRE%3d) %12.8f%12.8f%12.8f  CHARGE = %4.1f" \
     %(name, i+1, i+1, atom.coord[0]*f, atom.coord[1]*f, atom.coord[2]*f, atom.charge),

  # write basis
  basis = res.uncontracted_basis
  center = []
  type = []
  expo = []
  for prim in basis:
     c = prim.center
     for i,atom in enumerate(geom):
         atom_c = atom.coord
         if atom_c == c:
            atom_id = i+1
     sym = basis_correspond[prim.sym]
     center.append(atom_id)
     type.append(sym)
     expo.append(prim.expo)
  k=0
  for c in center:
    if k%20 == 0:
      print "\nCENTRE ASSIGNMENTS  ",
    print "%2d"%(c),
    k+=1
  k=0
  for t in type:
    if k%20 == 0:
      print "\nTYPE ASSIGNMENTS    ",
    print "%2d"%(t),
    k+=1
  k=0
  for e in expo:
    if k%5 == 0:
      print "\nEXPONENTS ",
    print ("%13.7e"%(e)).replace('e','D'),
    k+=1

  print '\n',

  # MOs
  if MO_type is None:
    MO_type = res.determinants_mo_type
  try:
    allMOs = res.uncontracted_mo_sets[MO_type]
  except KeyError:
    print "MO type not present in "+res.filename
    return
  mos = []
  occ = res.occ_num[MO_type]
  for i,mo in enumerate(allMOs):
    mos.append( (mo.eigenvalue, mo, i) )
  mos.sort()
  for i,orb2 in enumerate(mos):
   if occ[orb2[2]] > 0.:
    orb = orb2[1]
    print "MO %4d     MO 0.0        OCC NO =    %9.7f  ORB. ENERGY =%12.6f" \
      %(i+1, occ[orb2[2]], orb.eigenvalue),
    k=0
    for c in orb.vector:
      if k%5 == 0:
        print "\n",
      print "%15.8e"%(c),
      k+=1
    print "\n",
  print "END DATA"
  idx = res.mo_types.index(MO_type)
  try:
    print " THE  HF ENERGY =%20.12f THE VIRIAL(-V/T)="%(res.energies[idx]),
    print "%12.8f"%(res.virials[idx]) 
  except TypeError:
    print " THE  HF ENERGY =%20.12f THE VIRIAL(-V/T)=%12.8f"%(0.,0.)






def write_wfnFile(file,out=sys.stdout):
  if "basis" in file.mo_sets.keys() : 
    MoType = "basis"
  else:
    MoType = None
  file.clean_uncontractions()
  wfn_write(file,out,MoType)

  
fileTypes.append(wfnFile)

if __name__ == '__main__':
   main(wfnFile)

