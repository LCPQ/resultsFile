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



from math import *
a0 = 0.529177249

def prettyPrint(obj,shift=" "):
   if isinstance(obj,list):
     for k in obj:
       prettyPrint(k)
   elif isinstance(obj,dict):
     for k in obj:
       print k
       prettyPrint(obj[k])
   else:
     print str(obj)

def xproperty(fget, fset=None, fdel=None, doc=None):
   """Automatically defines the properties if the functions don't exist."""

   if isinstance(fget,str):
      attr_name = fget
      if fset is None:
         fset = fget
      fget = lambda obj: getattr(obj, attr_name)

   if isinstance(fset,str):
      attr_name = fset
      fset = lambda obj, val: setattr(obj, attr_name, val)

   return property(fget, fset, fdel, doc)



def build_get_funcs(var):
   line = ""
   line += """def get_"""+var+"""(self):
     try:
       getattr(self,'_"""+var+"""')
     except AttributeError:
       self._"""+var+""" = None
     return getattr(self, '_"""+var+"""')\n"""
   return line

def build_property(var,txt):
   line = "def _get_"+var+"(self):\n"
   line += "  return self.get_"+var+"()\n"
   line += var+" = property(_get_"+var+", fset=None, fdel=None, doc='"+txt+"')\n"
   return line

def get_data(var,txt,pos,type=""):
   line = """def get_%(var)s(self):
   if self._%(var)s is None:
      try:
        self.find_string("%(txt)s")
      except IndexError:
        return None
      pos = self._pos
      if pos is not None:
         line = self.text[pos].split()"""
   if type is "":
     line += """
         self._%(var)s = ' '.join(line[%(pos)s])"""
   else:
     line += """
         self._%(var)s = """+type+"""(' '.join(line[%(pos)s]))"""
   line += """
   return self._%(var)s"""
   line = line % {'var': var, 'txt': txt, 'pos': pos}
   return line

def init_variable(var,value):
  assert ( value is not None )
  line = """def get_"""+var+"""(self):
     try:
       getattr(self,'_"""+var+"""')
     except AttributeError:
       self._"""+var+""" = None
     if self._"""+var+""" is None:
       self._"""+var+""" = """+str(value)+"""
     return self._"""+var
  return line

def cartesianToSpherical(cart):
  x, y, z = cart
  r = sqrt(x**2+y**2+z**2)
  sign = 1.
  if x != 0.:
    sign = abs(x)/x
    theta = atan(y/x)
  else:
    sign = 1.
    theta = atan(pi/2.)
  if r != 0.:
    phi = acos(z/r)*sign
  else:
    phi = 0.
  return [r,theta,phi]

def sphericalToCartesian(sphe):
  r, theta, phi = sphe
  x = r*cos(theta)*sin(phi)
  y = r*sin(theta)*sin(phi)
  z = r*cos(phi)
  return [x,y,z]
  
  
