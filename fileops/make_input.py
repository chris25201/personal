#!/usr/bin/env python
#coding=utf-8

import os
import sys
import math


#for i in range(1, len(sys.argv)):
#  print i, sys.argv[i]


# how many disp POSCAR
s=sys.argv[1] 
ndisp=int(s)


for i in range(1,ndisp+1):
  s="%03d"%i
  poscar="POSCAR_"+s
  print(poscar)

  disp="disp"+s
  print(disp)
  if os.path.isdir(disp):
    print(disp+" already exist, skip it!")
    continue

  s="mkdir " + disp 
  os.system(s) 

  s="cp -f INCAR POTCAR KPOINTS vasp.lsf vdw_kernel.bindat " + disp
  os.system(s) 

  s="cp " + poscar+ " "  + disp + "/POSCAR"
  os.system(s) 

  s="cd "+disp + " ; bsub < vasp.lsf ; cd ../"
  os.system(s) 


