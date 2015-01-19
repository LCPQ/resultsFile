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



import struct
from library import *
import sys

class fortranBinary(object):
   """
Class for a fortran binary file.
Behaves like an array of records.
   """

   def __init__(self,name,mode):
       self._name = name
       self._mode = mode
       self._file = None
       self._recl = None
       self._form = None

   def get_file(self):
       if self._file is None:
          self._file = open(self._name,self._mode)
       return self._file

   def seek(self,index):
       self.file.seek(index)

   def next(self):
       if self.form is None:
          raise TypeError
       data = self.file.read(4)  # Rec length at the beginning of record
       if not data: raise IndexError
       self.recl = struct.unpack('l',data)[0]
       data = self.file.read(self._recl)
       data = struct.unpack(self.form,data)
       self.file.read(4)
       return list(data)

   def __iter__(self):
       if self.form is None:
          raise TypeError
       self.file.seek(0)
       while True:
          data = self.file.read(4)  # Rec length at the beginning of record
          if not data: raise StopIteration
          self.recl = struct.unpack('l',data)[0] 
          data = self.file.read(self._recl)
          data = struct.unpack(self.form,data)
          self.file.read(4)
          yield list(data)

   def __getitem__(self,index):
       if self.form is None:
          raise TypeError
       self.file.seek(0)
       for i in range(index-1):
          data = self.file.read(4)  # Rec length at the beginning of record
          if not data: raise IndexError
          recl = struct.unpack('l',data)[0] 
          data = self.file.read(recl+4)
       data = self.file.read(4)  # Rec length at the beginning of record
       if not data: raise IndexError
       self.recl = struct.unpack('l',data)[0] 
       data = self.file.read(self._recl)
       data = struct.unpack(self.form,data)
       self.file.read(4)
       return list(data)

   def close(self):
       self._file.close()


#  recl = xproperty ( '_recl', doc="Record length")
#  mode = xproperty ( '_mode', doc="R/W mode")
#  name = xproperty ( '_name', doc="R/W mode")
#  file = property  ( get_file, doc="Binary file")
#  form = xproperty ( '_form', doc="Format to read the binary record")
   for i in "recl mode name file form".split():
     exec """
def get_%(i)s(self): return self._%(i)s
def set_%(i)s(self,value): self._%(i)s = value
%(i)s = property(fget=get_%(i)s,fset=set_%(i)s) """%locals()



if __name__ == '__main__':
   f = fortranBinary(sys.argv[1],"rb")
   f.form = 'l'+4*15000*'h'+15000*'d'
   print f.name
   print f.mode
   print f[0]
   print f[10]
   print f[0]
