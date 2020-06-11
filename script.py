#!/usr/bin/env python
import os, sys, signal
from os import listdir
from os.path import isfile, join
from subprocess import call
from numpy import *
from random import shuffle

instType = "M";
myPath = "../Python/InstanceGenerator/" + instType + "/";

allFiles = [ f for f in listdir(myPath) if isfile(join(myPath,f)) ];
allFiles.sort();

for instance in allFiles:
    instanceName = instance.split(".")[0];
    noWaypoints = int(instanceName[7:9]);
    print "Running instance: " + instanceName;
    print "Number of waypoints:", noWaypoints;
    command = "./exeGTO " + myPath + instance + " 50 2";
    print command;
    var = os.system(command);
    if (var == 2):
        print "\nKeyboard interruption!\n";
        sys.exit();





