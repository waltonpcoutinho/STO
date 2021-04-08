#!/usr/bin/env python3
import os, sys, signal
from os import listdir
from os.path import isfile, join
from subprocess import call
from numpy import *
from random import shuffle

instanceType = str(sys.argv[1])

maxFleetSize = 2

myPath = "../Python/InstanceGenerator/" + instanceType + "/"
#myPath = "../GRTOP_instances/" + instanceType + "/"
allFiles = [ f for f in listdir(myPath) if isfile(join(myPath,f)) ]
allFiles.sort()

for instance in allFiles:
   instanceName = instance.split(".")[0]
   noWaypoints = int(instanceName[7:9])
   print("Running instance: " + instanceName)
   print("Number of waypoints:", noWaypoints)
   command = "./exeGTO " + myPath + instance + " 50 " + str(maxFleetSize) + "  " + instanceType + " STO "
   print(command)
   var = os.system(command);
   if (var == 2):
      print("\nKeyboard interruption!\n")
      sys.exit();





