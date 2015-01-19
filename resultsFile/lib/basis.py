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



#import pdb
from library import * 
from math import *
import sys

def normalize_basis_name(s):
  l = list(s)
  is_sphe = False
  l.sort()
  for i in range(10):
    if str(i) in l:
      is_sphe = True
  if is_sphe:
    if l[0] == '0':
      l.insert(0,'+')
    l = [l[-1]]+l[:-1]
  result = ""
  for i in l:
   result += i
  result = result.lower()
  return result


class primitive(object):
   """Class for a primitive basis function."""

   def __init__(self):
       _center = None  
       _expo   = None  
       _sym    = None  

   def __eq__(self,other):
       thr = 1.e-6
       result = True
       result = result and abs(self.center[0] - other.center[0]) < thr
       result = result and abs(self.center[1] - other.center[1]) < thr
       result = result and abs(self.center[2] - other.center[2]) < thr
       result = result and abs(self.expo - other.expo) < thr
       result = result and self.sym == other.sym
       return result 

   def __repr__(self):
     out = "%6s   %10.6f  %10.6f  %10.6f   %16.6e"%tuple(
        [self.sym]+list(self.center)+[self.expo])
     return out

   def __repr__debug__(self):
     out = ""
     out += str(self.sym)+'\t'
     out += str(self.center)+'\t'
     out += str(self.expo)
     return out

   def __cmp__(self,other):
       assert ( isinstance(other,primitive) )
       """Primitive functions are sorted according to the exponents."""
       if not isinstance(other,primitive):
          raise TypeError
       if self.expo < other.expo:
          return -1
       elif self.expo > other.expo:
          return 1
       elif self.expo == other.expo:
          return 0

   def get_norm(self):
       result = sqrt(self.overlap(self))
       return result

   norm = property (get_norm,None,doc="Sqrt( Integral f^2(R) dR ).")
   for i in "center expo sym".split():
     exec """
def get_%(i)s(self): return self._%(i)s
def set_%(i)s(self,value): self._%(i)s = value
%(i)s = property(fget=get_%(i)s,fset=set_%(i)s) """%locals()


#-----------------------

powersave = { 's':(0,0,0) }

def powers(sym):
  if sym in powersave:
     return powersave[sym]
  result = (sym.count('x'),sym.count('y'),sym.count('z'))
  powersave[sym] = result
  return result

  
fact_ = [1.]
def fact(n):
  global fact_
  nstart = len(fact_)
  if n < nstart :
    return fact_[n]
  else:
    for i in range(nstart,n+1):
      fact_.append(fact_[i-1]*float(i))
    return fact_[n]

def binom(n,m):
  return fact(n)/(fact(m)*fact(n-m))

def rintgauss(n):

  def ddfact2(n):
    if n%2 == 0: print 'error in ddfact2'
    res=1.
    for i in range(1,n+1,2):
      res*=float(i)
    return res

  res = sqrt(pi)
  if n == 0: return res
  elif n == 1: return 0.
  elif n%2 == 1: return 0.
  res /= 2.**(n/2)
  res *= ddfact2(n-1)
  return res

          
def Goverlap(fA,fB):
  fA = spheToCart(fA)
  fB = spheToCart(fB)
  if isinstance(fA,contraction):
    result = fA.overlap(fB)
  elif isinstance(fB,contraction):
    result = fB.overlap(fA)
  else:
    result = GoverlapCart(fA,fB)
  return result
          

def GoverlapCart(fA,fB):
  gamA=fA.expo
  gamB=fB.expo
  gamtot = gamA+gamB
  SAB=1.0
  A = fA.center
  B = fB.center
  nA = powers(fA.sym)
  nB = powers(fB.sym)
  for l in range(3):
    Al = A[l]
    Bl = B[l]
    nAl = nA[l]
    nBl = nB[l]
    u=gamA/gamtot*Al+gamB/gamtot*Bl
    arg=gamtot*u*u-gamA*Al*Al-gamB*Bl*Bl
    alpha=exp(arg)/gamtot**((1.+float(nAl)+float(nBl))/2.)
    temp = sqrt(gamtot)
    wA=temp*(u-Al)
    wB=temp*(u-Bl)
    accu=0.
    for n in range (nAl+1):
      wAn = wA**n * binom(nAl,n)
      for m in range (nBl+1):
        integ=nAl+nBl-n-m
        accu+=wAn*wB**m*binom(nBl,m)*rintgauss(integ)
    SAB*=accu*alpha
  return SAB

def GoverlapCartNorm2(fA,fB):
  gamA=fA.expo
  gamB=fB.expo
  gamtot = gamA+gamB
  SAB=1.0
  nA = powers(fA.sym)
  nB = powers(fB.sym)
  for l in range(3):
    nAl = nA[l]
    nBl = nB[l]
    SAB*=rintgauss(nAl+nBl)/(gamA+gamB)**((1.+float(nAl)+float(nBl))/2.)
  return SAB

def GoverlapCartNorm(fA):
  gamA=fA.expo
  SAB=1.0
  nA = powers(fA.sym)
  for l in range(3):
    nAl = nA[l]
    SAB*=rintgauss(2*nAl)/(2.*gamA)**(0.5+float(nAl))
  return SAB

angular_momentum = {}
for i,l in enumerate("spdfghijklmno"):
  angular_momentum[l]=i

def get_lm(sym):
  if 'x' in sym or 'y' in sym or 'z' in sym or sym == 's':
    return None, None
  else:
    return angular_momentum[sym[0]], int(sym[1:])

def xyz_from_lm(l,m):
  """Returns the xyz powers and the coefficients of a spherical function
     expressed in the cartesian basis"""
  power = []
  coef  = []
  absm = abs(m)
  nb2  = absm
  nb1  = (l-absm)/2
  clmt = [ (-0.25)**t * binom(l,t) * binom(l-t,absm+t) for t in xrange(nb1+1) ]
  mod_absm_2 = absm % 2
  if m>=0:
    nb2_start = mod_absm_2
  elif m<0:
    nb2_start = (absm+1) % 2
  if m != 0:
    norm = 2.**(-absm)/fact(l) * sqrt(2.*abs(fact(l+m)*fact(l-m)))
  else:
    norm = 1.
  for n1 in range(nb2_start,nb2+1,2):
    k = (absm-n1)/2
    factor = (-1.)**k * binom(absm,n1) * norm
    for t in xrange(nb1+1):
      for n2 in xrange(t+1):
        coe = clmt[t] * factor * binom(t,n2)
        ipw = ( n1+2*n2 , absm-n1+2*(t-n2), l-2*t-absm )
        done = False
        for i,c in enumerate(coef):
         if c != 0.:
          if not done and ipw in power:
            idx = power.index(ipw)
            coef[idx] += coe
            done = True
        if not done:
          power.append(ipw)
          coef.append(coe)
  return power,coef

def spheToCart(fA):
  l,m = get_lm(fA.sym)
  if l is None:
    return fA
  power, coef = xyz_from_lm(l,m)

  contr = contraction()
  for p,c in zip(power,coef):
    gauss = gaussian()
    gauss.center = fA.center
    gauss.expo = fA.expo
    gauss.sym  = ''
    for l,letter in enumerate('xyz'):
      gauss.sym += p[l]*letter
    contr.append(c,gauss)
  return contr
    
        
class gaussian(primitive):
   """Gaussian primitive function."""

   def __init__(self):
     primitive.__init__(self)

   def overlap(self,other):
     result = Goverlap(self,other)
     return result

   def value(self,r):
     """Value at r."""
     x, y, z = r
     x -= self.center[0]
     y -= self.center[1]
     z -= self.center[2]
     r2 = x**2 + y**2 + z**2
     px, py, pz = powers(self.sym)
     P = x**px * y**py * z**pz
     return P*exp(-self.expo*r2)

   def valuer2(self,r2):
     """Value at sqrt(r2)."""
     return exp(-self.expo*r2)

#-----------------------

class contraction(object):
   """Contraction of primitive functions."""

   def __init__(self):
       self._center = None
       self._prim = []
       self._coef = []
       self._sym  = None

   def __eq__(self,other):
     thr = 1.e-6
     result = True
     result = result and abs(self._center[0] - other._center[0]) < thr
     result = result and abs(self._center[1] - other._center[1]) < thr
     result = result and abs(self._center[2] - other._center[2]) < thr
     result = result and self._sym == other._sym
     result = result and len(self._coef) == len(other._coef)
     if result:
       for i in range(len(self._coef)):
         result = result and abs(self._coef[i] - other._coef[i]) < thr
     if result:
       for i in range(len(self._prim)):
         result = result and self._prim[i] == other._prim[i]
     return result
      
   def __repr__(self):
       out = "%6s   %10.6f  %10.6f  %10.6f\n"%tuple(
        [self.sym]+list(self.center))
       for i,a in enumerate(self.prim):
         out += "  %16.6e %16.6e\n"%(a.expo,self.coef[i])
       return out

   def append(self,c,func):
       if self.sym is None:
          self.sym = func.sym
       if self.center is None:
          self.center = func.center
       assert ( isinstance(func,primitive) )
       assert ( self.center == func.center )
       self._prim.append(func)
       self._coef.append(c)

   def get_norm(self):
       sum=0.
       for i, ci in enumerate(self.coef):
         for j, cj in enumerate(self.coef):
           sum += ci*cj*self.prim[i].overlap(self.prim[j])/(self.prim[i].norm*self.prim[j].norm)
       result = sqrt(sum)
       return result

   def overlap(self,other0):
       if isinstance(other0,contraction):
         other=other0
       else:
         other=contraction()
         other.sym = other0.sym
         other.center = other0.center
         other.append(1.,other0)
       prim = self.prim
       oprim = other.prim
       sum = 0.
       for i, ci in enumerate(self.coef):
         for j, cj in enumerate(other.coef):
           sum += ci*cj*prim[i].overlap(oprim[j])
       return sum

   def get_coef(self):
     return self._coef
   coef = property(fget=get_coef)

   def get_center(self):
     return self._center
   center = property(fget=get_center)

   def get_prim(self):
     return self._prim
   prim = property(fget=get_prim)
   def normalize(self):
       coef = self.coef
       prim = self.prim
       n = self.norm
       for i, ci in enumerate(coef):
           coef[i] /= n

   def value(self,r):
     """Value at r."""
     x, y, z = r
     x -= self.center[0]
     y -= self.center[1]
     z -= self.center[2]
     r2 = x**2 + y**2 + z**2
     px, py, pz = powers(self.sym)
     P = x**px * y**py * z**pz
     result = 0.
     for c,chi in zip(self.coef,self.prim):
       result += c*chi.valuer2(r2)
     return P*result

            
   for i in "center prim sym coef".split():
     exec """
def get_%(i)s(self): return self._%(i)s
def set_%(i)s(self,value): self._%(i)s = value
%(i)s = property(fget=get_%(i)s,fset=set_%(i)s) """%locals()

   norm = property (get_norm,None,doc="Sqrt( Integral f^2(R) dR ).")


try:
  import lib_cython
  fact = lib_cython.fact
  binom = lib_cython.binom
  rintgauss = lib_cython.rintgauss
  powers = lib_cython.powers
  GoverlapCart = lib_cython.GoverlapCart
  GoverlapCartNorm2 = lib_cython.GoverlapCartNorm2
except ImportError:
  pass

#--------------

if __name__ == '__main__':
 for l in range(6):
   print '---'
   for m in range(-l,l+1):
     print "%3d %3d :"%(l,m), xyz_from_lm(l,m)
 sys.exit(0)

if __name__ == '__main__':
        a = gaussian()
        a.sym = 'f+0'
        a.expo = 3.
        a.center = [0.,0.,0.]
        b = gaussian()
        b.sym = 'g+2'
        b.expo =  4.
        b.center = [0.,0.,0.]
        x = spheToCart(a)
        y = spheToCart(a)
        print x
        sys.exit(0)
        print GoverlapCart(a,b)
        print GoverlapCartNorm2(a,b)
        print ''
        for i in range(0,10):
          print rintgauss(i)
        sys.exit(0)
        c = gaussian()
        c.sym = 's'
        c.expo = 1.159
        c.center = [0.,0.,0.]
        d = contraction()
        d.append(0.025494863235, a)
        d.append(0.190362765893, b)
        d.append(0.852162022245, c)
        d.normalize()
        print d
        sys.exit(0)
        a = gaussian()
        a.sym = 'x'
        a.expo = 1.0
        a.center = [0.,0.,0.]
        b = gaussian()
        b.sym = 'x'
        b.expo = 2.0
        b.center = [0.,0.,0.]
        c = contraction()
        c.append(1.,a)
        c.append(2.,b)
        print c.value( (1.,1.,0.) )
        print b.overlap(a)
        print a.overlap(b)
        print a.norm
        print c.norm, c.coef
        c.normalize()
        print c.norm, c.coef
        c.normalize()
        print c.norm, c.coef
        print Goverlap(a,b)
 
