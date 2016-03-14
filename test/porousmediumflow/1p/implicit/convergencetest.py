#!/usr/bin/env python

from math import *
import subprocess
import sys

if len(sys.argv) != 2:
    #sys.cerr.write('Please provide a single argument <testname> to the script')
    print('Please provide a single argument <testname> to the script')

testname = str(sys.argv[1])

subprocess.call(['rm', testname + '.log'])
print("Removed old log file ({})!".format(testname + '.log'))
for i in [10, 100, 1000, 5000, 10000, 15000, int(1e5)]:
    subprocess.call(['./' + testname,
                     '-Grid.Cells', str(i)+" 1",
                     '-Problem.Name', testname])

logfile = open(testname + '.log', "r+")

error = []
hmax = []
for line in logfile:
    line = line.strip("\n")
    line = line.strip("\[ConvergenceTest\]")
    line = line.split()
    error.append(float(line[2]))
    hmax.append(float(line[5]))

logfile.truncate(0)
for i in range(len(error)-1):
    if not (error[i] < 1e-12 or error[i] < 1e-12):
        rate = (log(error[i])-log(error[i+1]))/(log(hmax[i])-log(hmax[i+1]))
        message = "Error: " + str(error[i]) + " hMax: " + str(hmax[i]) + " Rate: " + str(rate) + "\n"
        logfile.write(message)
    else:
        message = "error = 0 (exact solution)"
        logfile.write(message)

logfile.close()
print("\nComputed the following convergence rates:")
subprocess.call(['cat', testname + '.log'])
