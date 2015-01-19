#!/usr/bin/python
from math import *

powersave = { 's':(0,0,0) }
def powers(sym):
  if sym in powersave:
     return powersave[sym]
  result = (sym.count('x'),sym.count('y'),sym.count('z'))
  powersave[sym] = result
  return result

fact_ = [1.]
cpdef fact(int n):
  global fact_
  cdef int nstart = len(fact_)
  cdef int i
  if n >= nstart :
    for i in range(nstart,n+1):
      fact_.append(float(fact_[i-1]*i))
  return fact_[n]

def binom(int n,int m):
  return fact(n)/(fact(m)*fact(n-m))


cdef ddfact2(int n):
  if n%2 == 0: print 'error in ddfact2'
  cdef double res=1.
  cdef int i
  for i in range(1,n+1,2):
    res*=float(i)
  return res

cdef double sqpi = sqrt(pi)

cpdef rintgauss(int n):
  res = sqpi
  if n == 0: return res
  elif n == 1: return 0.
  elif n%2 == 1: return 0.
  res /= 2.**(n/2)
  res *= ddfact2(n-1)
  return res

cpdef GoverlapCart(fA,fB):
  cdef double gamA=fA.expo
  cdef double gamB=fB.expo
  cdef double gamtot = gamA+gamB
  cdef double SAB=1.0
  cdef int l, n, m
  cdef double u, arg, alpha, temp, wA, wB, accu
  cdef int integ
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
    alpha=exp(arg)/gamtot**((1.+float(nAl)+float(nBl))*0.5)
    temp = sqrt(gamtot)
    wA=temp*(u-Al)
    wB=temp*(u-Bl)
    accu=0.
    for n in range (nAl+1):
      for m in range (nBl+1):
        integ=nAl+nBl-n-m
        accu+=wA**n*wB**m*binom(nAl,n)*binom(nBl,m)*rintgauss(integ)
    SAB*=accu*alpha
  return SAB

cpdef GoverlapCartNorm2(fA,fB):
  cdef double gamA=fA.expo
  cdef double gamB=fB.expo
  cdef double gamtot = gamA+gamB
  cdef double SAB=1.0
  cdef int l
  nA = powers(fA.sym)
  nB = powers(fB.sym)
  for l in range(3):
    nAl = nA[l]
    nBl = nB[l]
    SAB*=rintgauss(nAl+nBl)/(gamA+gamB)**((1.+float(nAl)+float(nBl))/2.)
  return SAB

