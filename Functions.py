from ConfigParser import SafeConfigParser
import sys, os, traceback, time
import subprocess as sp
from os.path import isfile, join
# heck if a folder exists
def CheckFolderExist(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)
# Function gets the list of files with the extension FileExtension, in the specific folder PathFile
def GetListFile(PathFile, FileExtension):
    return [os.path.splitext(f)[0] for f in os.listdir(PathFile) if
            isfile(join(PathFile, f)) and os.path.splitext(f)[1] == '.' + FileExtension]
#function to cretae a log file
class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
def timeStamp():
    t = time.localtime()
    month = t.tm_mon
    day = t.tm_mday
    hour = t.tm_hour
    minute = t.tm_min
    second = t.tm_sec
    return "at %i:%i:%i, %i/%i"%(hour, minute, second, month, day)


def spawnProcesses(maxProcesses, processNames, argListDict):
# spawn a given number of child processes (spawning more as children complete)
# and return stdouts and stderrs
    processNamesLocal = list(processNames)
    #[os.path.basename(name) for name in list(processNames)]
    processes = {}
    stdOuts = {}
    stdErrs = {}
    stdOutFiles = {}
    stdErrFiles = {}
    while len(processNamesLocal) != 0 or len(processes) != 0:
        time.sleep(1) # make this a less-busy wait
        if len(processes) < maxProcesses and len(processNamesLocal) > 0:
            procName = processNamesLocal.pop()
            #print "procName = "+procName
            #print "argList = "+str(argListDict[procName])
            stdOutFiles[procName] = open(os.path.join(tempDir,os.path.basename(procName)+".stdout"),"w")
            #print "Created temp file %s for writing."%(os.path.join(tempDir,os.path.basename(procName)+".stdout"))
            stdErrFiles[procName] = open(os.path.join(tempDir,os.path.basename(procName)+".stderr"),"w")
            # TODO: add timestamp to each stdout as it completes, instead of logging time that all processes complete?
            # or maybe define an output handling callback function and call it after each process completes
            processes[procName] = sp.Popen(argListDict[procName], stdout=stdOutFiles[procName], stderr=stdErrFiles[procName])
        if len(processes.keys()) != 0:
            for procName in processes.keys():
                if processes[procName].poll() is not None:
                    processes.pop(procName)
    for procName in stdOutFiles.keys():
        stdOutFiles[procName].close()
        stdErrFiles[procName].close()
        stdOutFile = open(os.path.join(tempDir,os.path.basename(procName)+".stdout"),"rU")
        #print "Opened temp file %s for reading."%(os.path.join(tempDir,os.path.basename(procName)+".stdout"))
        stdErrFile = open(os.path.join(tempDir,os.path.basename(procName)+".stderr"),"rU")
        stdOuts[procName] = stdOutFile.read()
        stdErrs[procName] = stdErrFile.read()
        stdOutFile.close()
        stdErrFile.close()
    return stdOuts, stdErrs
tempDir = os.path.abspath(os.path.join(os.getcwd(),"temp"))
CheckFolderExist(tempDir)