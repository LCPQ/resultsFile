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



"""resultsFile library."""

__author__ = "Anthony SCEMAMA <scemama@irsamc.ups-tlse.fr>"
__date__   = "25 Oct 2007"

import os

wd = os.path.dirname(__file__)
all = [ i[:-3] for i in os.listdir(wd) if i.endswith("File.py") ]
if "q5CostFile" in all:
  all.remove("q5CostFile")
  all += [ "q5CostFile" ]

for mod in all:
   try:
     exec 'from '+mod+' import *'
   except:
     print "Error importing module", mod
     pass


