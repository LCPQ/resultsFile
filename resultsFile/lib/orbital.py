#!/usr/bin/env python 
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



from library import * 
from math import *

class orbital(object):
   """Class for an orbital."""

   def __init__(self):
       self._set = None
       self._eigenvalue = None
       self._vector = []
       self._basis = None
       self._sym = None
       
   def __repr__(self):
       out = "%10s %6s %10.6f\n"%(\
          self.set, self.sym, self.eigenvalue)
       count=0
       for c in self.vector:
         out += "%17.8e "%(c,)
         count += 1
         if count == 5:
            count = 0
            out += '\n'
       return out

   def __repr__debug__(self):
       out = ""
       out += "Orbital:\n"
       out += " Set         : "+str(self.set)+'\n'
       out += " Sym         : "+str(self.sym)+'\n'
       out += " Basis       : "+str(self.basis)+'\n'
       out += " Eigenvalue  : "+str(self.eigenvalue)+'\n'
       out += " Vector      : "+str(self.vector)+'\n'
       return out

   def __cmp__(self,other):
       assert ( isinstance(other,orbital) )
       if self.eigenvalue < other.eigenvalue:
          return -1
       elif self.eigenvalue > other.eigenvalue:
          return 1
       elif self.eigenvalue == other.eigenvalue:
          return 0

   def dot(self,other):
       assert ( isinstance(other,orbital) )
       result = 0.
       for i,j in zip(self.vector,other.vector):
           result += i*j
       return sqrt(result)

   def is_ortho(self,other):
       assert ( isinstance(other,orbital) )
       if self.dot(other) == 0. :
          return True
       else:
          return False

   def norm(self):
       return self.dot(self)

   def value(self,r):
     """Value at r."""
     x, y, z = r
     result = 0.
     for c,chi in zip(self.vector,self.basis):
       result += c*chi.value(r)
     return result

   for i in "eigenvalue vector basis set sym".split():
     exec """
def get_%(i)s(self): return self._%(i)s
def set_%(i)s(self,value): self._%(i)s = value
%(i)s = property(fget=get_%(i)s,fset=set_%(i)s) """%locals()





if __name__ == '__main__':
   l = []
   l.append(orbital())
   l.append(orbital())
   l.append(orbital())
   l[0].eigenvalue = 2.
   l[1].eigenvalue = 1.
   l[2].eigenvalue = 3.
   for i in l:
     print i.eigenvalue
   l.sort()
   print ""
   for i in l:
     print i.eigenvalue
   l[0].vector = [1., 0., 0.]
   l[1].vector = [0., 1., 0.]
   l[2].vector = [0., 1., 1.]
   print l[0].is_ortho(l[1])
   print l[1].is_ortho(l[2])
   print l[1].dot(l[2])
   print l[1].dot(l[0])
   print l[1].norm()
   print l[2].norm()



