import os
import sys
import math


#for i in range(1, len(sys.argv)):
#  print i, sys.argv[i]


# how many disp POSCAR
s=sys.argv[1]
ncal=int(s)


for i in range(1,ncal+1):
  s="%03d"%i

  cal="cal"+s
  print(cal)
  if os.path.isdir(cal):
    print(cal+" already exist, skip it!")
    continue

  s="mkdir " + cal
  os.system(s)

  s="cp -f nep.txt run.in gpumd.sh model.xyz " + cal
  os.system(s)

  s="cd "+ cal + " ; python ../randomremoveO.py ; bsub < gpumd.sh ; cd ../"
  os.system(s)


