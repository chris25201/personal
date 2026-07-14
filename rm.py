#!/usr/bin/env python
#coding=utf-8

import os
import sys
import math


#for i in range(1, len(sys.argv)):
#  print i, sys.argv[i]


# how many disp POSCAR
n=sys.argv[1]
m=str(int(n))
con=["222","331","332","442"]
os.system("cp -r ../v24/rm" + n + " rm" + n)
os.chdir("rm" + n)
os.system("pwd; rm *xyz")

for i in range(0,4):
  f=con[i]

#  if os.path.isdir(disp):
#    print(disp+" already exist, skip it!")
#    continue

  os.chdir(f)
  os.system("pwd")

  s="rm -r disp* POSCAR_*  struct* movie* nep.txt NEP* err* out* log* *out make_dataset.xyz"
  os.system(s)
#  print(s)
#  os.system("pwd")

  s="cp ../../nep.txt . "
  os.system(s)

  os.chdir("..")
#  s="cd ../../../v" + n +"/md_v" + n + "/" + f +  "; bsub < vasp.lsf ; cd ../.."
#  os.system(s)


