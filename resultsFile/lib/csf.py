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

class CSF(object):
   """Class for an configuration state function."""

   def __init__(self):
       self._determinants = []
       self._coefficients = []
       
   def __repr__(self):
       out = ""
       for d,c in zip(self.determinants,self.coefficients):
         out += "%10.8f\n"%(c,)
         for spin in d:
           for orb in d[spin]:
             out += "%4d "%(orb,)
           out += "\n"
       return out

   def __repr__debug__(self):
       out = ""
       out += "CSF:\n"
       out += " Determinants: "+str(self.determinants)+'\n'
       out += " Coefficients: "+str(self.coefficients)+'\n'
       return out

   def __cmp__(self,other):
       assert ( isinstance(other,CSF) )
       if len(self.coefficients) < len(other.coefficients):
          return -1
       elif len(self.coefficients) > len(other.coefficients):
          return 1
       elif len(self.coefficients) == len(other.coefficients):
          return 0

   def append(self,c,up,dn):
       det = {}
       det['alpha'] = up
       det['beta'] = dn
       self._determinants.append( det )
       self._coefficients.append( c )

   for i in "determinants coefficients".split():
     exec """
def get_%(i)s(self): return self._%(i)s
def set_%(i)s(self,value): self._%(i)s = value
%(i)s = property(fget=get_%(i)s,fset=set_%(i)s) """%locals()

