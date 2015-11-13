#!/usr/bin/python
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


n = range(5,35)
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

plt.plot(n,mass)
plt.xlabel('number of nodes')
if sys.argv[2]=='time':
  plt.ylabel('manouvre time (s)')
else:
  plt.ylabel('propellant mass (kg)')

plt.figure()
plt.plot(n,cputime,'o')
plt.xlabel('number of nodes')
plt.ylabel('CPU time (ms)')
plt.xlim([6,50])
plt.ylim([0,1000])
plt.show()
