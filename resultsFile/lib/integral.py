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



from library import * 
from math import *

class integral(object):
   """Class for an integral."""

   def __init__(self):
       self._indices = None
       self._value = None
       
   def __repr__(self):
       out = ""
       i = self.indices
       if len(i) == 4:
         out += "( %4d %4d | %4d %4d )"% tuple( [int(i[k]) for k in range(4)] )
       elif len(i) == 2:
         out += "( %4d | h | %4d )"% tuple( [int(i[k]) for k in range(2)] )
       out += " : %24.16e" % ( self.value )
       return out

   def __cmp__(self,other):
       assert ( isinstance(other,integral) )
       if self._value < other._value:
          return -1
       elif self._value > other._value:
          return 1
       elif self._value == other._value:
          return 0

   for i in "indices value".split():
     exec """
def get_%(i)s(self): return self._%(i)s
def set_%(i)s(self,value): self._%(i)s = value
%(i)s = property(fget=get_%(i)s,fset=set_%(i)s) """%locals()


if __name__ == '__main__':
   i = integral()
   i.indices = [1,4]
   i.value = 1.5
   print i
   j = integral()
   j.indices = [4,3,2,1]
   j.value = 2.5
   print j
   print i<j
