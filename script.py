#!/usr/bin/env python3
import os, sys, signal
from os import listdir
from os.path import isfile, join
from subprocess import call
from numpy import *
from random import shuffle

#get instance type from command line
instanceType = str(sys.argv[1])

#remove results file 
myres = "Results/compResults_" + instanceType + ".out"
if os.path.exists(myres):
   os.remove(myres)

#remove error log
myerr = "Results/runLog_" + instanceType + ".out"
if os.path.exists(myerr):
   os.remove(myerr)

myPath = "../Python/InstanceGenerator/" + instanceType + "/"
allFiles = [ f for f in listdir(myPath) if isfile(join(myPath,f)) ]
allFiles.sort()

for instance in allFiles:
   instanceName = instance.split(".")[0]
   noWaypoints = int(instanceName[7:9])
   maxFleetSize = floor(noWaypoints/2)
   print("Running instance: " + instanceName)
   print("Number of waypoints:", noWaypoints)
   command = "./exeGTO " + myPath + instance + " 50 " + str(maxFleetSize) + " 2>&1 | tee -a Results/runLog_" + instanceType + ".out"
   print(command)
   var = os.system(command);
   if (var == 2):
      print("\nKeyboard interruption!\n")
      sys.exit();





