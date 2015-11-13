#!/usr/bin/python
#Tests a model against different number of nodes
#EXAMPLE: python test main_simple.mod mass
import subprocess, sys
from matplotlib import pyplot as plt
import time

if len(sys.argv)==1:
  raise ValueError("Specify which model to use")

f = open(sys.argv[1],'r')
lines = f.readlines()
f.close()
mass = []
cputime = []


n = range(6,20)
for i in n:
  lines[13] = '\tparam n:=' + str(i) + ';\t\t\t\t\t#Numbers of nodes\r\n' 
  f = open(sys.argv[1],'w')
  f.writelines(lines)
  f.close()
  start = time.time()
  subprocess.call(["/opt/ampl/ampl",sys.argv[1]])
  end = time.time()
  f = open('out/sol.out','r')
  out = f.readlines()
  f.close()
  if sys.argv[2]=='time':
    mass.append(eval(out[-1].split()[0]))
  else:
    mass.append(9472.06 - eval(out[-1].split()[1]))
  cputime.append((end-start)*1000)

seg = range(n[0]-2,n[-1]-1)
plt.plot(seg,mass)
plt.xlabel('number of segments')
if sys.argv[2]=='time':
  plt.ylabel('manouvre time (s)')
else:
  plt.ylabel('propellant mass (kg)')

plt.figure()
plt.plot(seg,cputime,'o')
plt.xlabel('number of segments')
plt.ylabel('CPU time (ms)')
plt.xlim([6,50])
plt.ylim([0,1000])
plt.show()
