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


from lib import *
import lib.basis as Basis

def get_uncontracted_basis(basis):
   uncontr = []
   for contr in basis:
      for b in contr.prim:
         uncontr.append(b)
   return uncontr

def get_uncontracted_mo_sets(basis,uncontracted_basis,mo_sets,mo_types):
   cdef dict uncontr = {}
   cdef int lenmovector
   cdef int lenbasis = len(basis)
   cdef double ci
   cdef int i, j
   for motype in mo_types:
     uncontr[motype] = []
     for mo in mo_sets[motype]:
       lenmovector = len(mo.vector)
       monew = orbital()
       monew.basis = uncontracted_basis
       monew.eigenvalue = mo.eigenvalue
       monew.set = motype
       v = []
       for i in range(lenbasis):
         contr = basis[i]
         if i<lenmovector:
           ci = mo.vector[i]
           if ci == 0.:
             for j in range(len(contr.prim)):
                v.append(0.)
           else:
             for p, c in zip(contr.prim,contr.coef):
                v.append(c*ci/p.norm)
       monew.vector = v
       uncontr[motype].append(monew)
   return uncontr


def clean_contractions(basis):
  newbasis = []
  cdef int i, k, l, lenmovector
  idx = range(len(basis))
  for k,b1 in enumerate(basis):
    addBasis=True
    for l, b2 in enumerate(basis[:k]):
      if b2 == b1:
        idx[k] = l
        addBasis=False
        break
    if addBasis:
      newbasis.append(b1)
  self._basis = newbasis 

  mo_sets = self.mo_sets
  for motype in self.mo_types:
    for mo in mo_sets[motype]:
      lenmovector = len(mo.vector)
      newvec = [None for i in idx]
      for i in idx:
        newvec[i] = 0.
      for k,l in enumerate(idx):
        if k < lenmovector:
          newvec[l] += mo.vector[k]
      mo.vector = []
      for c in newvec:
        if c is not None:
          mo.vector.append(c)

#  def clean_uncontractions(self):
#    basis = self.uncontracted_basis
#    newbasis = []
#    idx = range(len(basis))
#    for k,b1 in enumerate(basis):
#      addBasis=True
#      for l, b2 in enumerate(basis[:k]):
#        if b2 == b1:
#          idx[k] = l
#          addBasis=False
#          break
#      if addBasis:
#        newbasis.append(b1)
#    self._uncontracted_basis = newbasis 

#    mo_sets = self.uncontracted_mo_sets
#    for motype in self.mo_types:
#      for mo in mo_sets[motype]:
#        lenmovector = len(mo.vector)
#        newvec = [None for i in idx]
#        for i in idx:
#          newvec[i] = 0.
#        for k,l in enumerate(idx):
#          if k < lenmovector:
#            newvec[l] += mo.vector[k]
#        mo.vector = []
#        for c in newvec:
#          if c is not None:
#            mo.vector.append(c)
#      
#  def convert_to_cartesian(self):
#    basis = self.basis
#    newbasis = []
#    idx = range(len(basis))
#    map = []
#    weight = []
#    for i,b in enumerate(basis):
#      l, m = Basis.get_lm(b.sym)
#      if l is None:
#        newbasis.append(b)
#        map.append(i)
#        weight.append(1.)
#      else:
#        powers, coefs = xyz_from_lm(l,m)
#        for j,prim in enumerate(b.prim):
#          b.coef[j] /= prim.norm
#        for c, p in zip(coefs, powers):
#          contr = copy.deepcopy(b)
#          sym = ''
#          for l,letter in enumerate('xyz'):
#            sym += p[l]*letter
#          contr.sym = sym
#          for j,prim in enumerate(contr.prim):
#            prim.sym = sym
#            contr.coef[j] *= prim.norm
#          newbasis.append(contr)
#          map.append(i)
#          weight.append(c)

#    mo_sets = self.mo_sets
#    for motype in self.mo_types:
#      for mo in mo_sets[motype]:
#        newvec = []
#        vec = mo.vector
#        for i,w in zip(map,weight):
#          newvec.append(vec[i]*w)
#        mo.vector = newvec

#    same_as = {}
#    for i,b1 in enumerate(newbasis):
#      for j,b2 in enumerate(newbasis[:i]):
#        if b1 == b2:
#         same_as[i] = j
#         weight[j] += weight[i]
#         break
#    to_remove = same_as.keys()
#    to_remove.sort()
#    to_remove.reverse()
#    for i in to_remove:
#      newbasis.pop(i)
#      weight.pop(i)
#      map.pop(i)


#    for motype in self.mo_types:
#      for mo in mo_sets[motype]:
#        for i in to_remove:
#          index = same_as[i]
#          value = mo.vector.pop(i)
#          mo.vector[index] += value

#    self._basis = newbasis
#    self._mo_sets = mo_sets

#  def find_string(self,chars):
#      """Finds the 1st occurence of chars.
#      """
#      self._pos = 0
#      self.find_next_string(chars)

#  def find_last_string(self,chars):
#      """Finds the 1st occurence of chars.
#      """
#      self._pos = len(self.text)-1
#      self.find_prev_string(chars)

#  def find_next_string(self,chars):
#      """Finds the next occurence of chars.
#      """
#      pos  = self._pos
#      text = self.text
#      found = False
#      while not found and pos < len(text):
#         if chars in text[pos]:
#            found = True
#         else:
#            pos += 1
#      if not found:
#         raise IndexError
#      self._pos = pos

#  def find_prev_string(self,chars):
#      """Finds the next occurence of chars.
#      """
#      pos  = self._pos
#      text = self.text
#      found = False
#      while not found and pos < len(text):
#         if chars in text[pos]:
#            found = True
#         else:
#            pos -= 1
#      if not found:
#         raise IndexError
#      self._pos = pos


#  for i, j in local_vars:
#     if i not in defined_vars:
#        exec build_get_funcs(i) in locals()
#     exec build_property(i,j) in locals()
#  del i,j


