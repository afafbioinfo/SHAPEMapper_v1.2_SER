#!/usr/bin/env python
from ConfigParser import SafeConfigParser
from os.path import isfile, join
import sys, os, traceback, time
import subprocess as sp
import Functions as Fct




config= SafeConfigParser()
config.read("Shape.cfg")
#Get the ensemble of folders path
outputDir=os.path.join( os.getcwd(),config.get("Paths","outputdir"))
fastaDir=os.path.join(outputDir,config.get("Paths","fastaDir"))
alignmentDir=os.path.join(outputDir,config.get("Paths","alignmentDir"))
removeAmbigDel=config.get("Parameters","removeAmbigDel")
minphredtocount=float(config.get("Parameters","minphredtocount"))
maxProc=config.get("Parameters","maxProc")
FastaFileExtension=config.get("Parameters","FastaFileExtension")
Experiments=(config.get("Parameters","Experiments")).split(',')

#Redirect all the print to Logfile.txt
sys.stdout =Fct.Logger("Logfile.txt")

#print Experiments[0],Experiments[1]
ListRNAs=Fct.GetListFile(fastaDir, FastaFileExtension)
#print ListRNAs
'''
#1 Run parseAlignment.cpp Parse aligned reads
try:
    alignedFileList = os.listdir(alignmentDir)
    samFilePathDict = {} # sample name:aligned file path
    indexPathDict = {} # sample name:reference fasta file path
    for filee in alignedFileList:
        sampleName=filee.split(".")[0]
        #print sampleName
        samFilePathDict[sampleName] = os.path.join(alignmentDir,filee)
        indexPathDict[sampleName] = os.path.join( fastaDir,sampleName+"."+FastaFileExtension)

        argListDict = {}
        options = ["-combine_strands","-deletion_masking","-randomly_primed","-trim_both_ends"]
        selected = ["True", "False", "off", "off"]
        argList = ["./parseAlignment"]
        for i in range(len(options)):
                if selected[i] == True:
                    argList.append(options[i])
        mutationstringDir=os.path.join(outputDir,"mutation_strings/")
        argList += ["-primer_length", str(10),
        "-min_map_qual", str(30),
         "-ref_seqs", indexPathDict[sampleName],
        "-file_in",samFilePathDict[sampleName],
         "-out_folder",mutationstringDir ]
        if removeAmbigDel==True:
                argList += ["-remove_ambig_del"]
        argListDict[sampleName] = argList

        sampleNameList=[sampleName]
        stdOuts, stdErrs =Fct.spawnProcesses(maxProc, sampleNameList, argListDict)
        #print"*",indexPathDict[sampleName]
        for sampleName in stdErrs:
                stdErr = stdErrs[sampleName]
                stdOut = stdOuts[sampleName]
                if stdOut != "" and stdOut is not None:
                    print str(stdOut)+"\n"

except:
    errorString = "Error:" + traceback.format_exc() + "\n"
    errorString += "Alignment parsing failed %s."%Fct.timeStamp()

    print errorString
    sys.exit(1)

#-------------------------------------------------------------------------------------------
#2 Run Count mutations

try:
    print "\nStarting mutation counting at %s\n"%(Fct.timeStamp())
    # get list of parsed mutation string files
    #mutationStringDir = os.path.join(outputDir,"mutation_strings/")
    folderFileList = os.listdir(mutationstringDir)
    argListDict = {}
    txtFileList = []
    for sampleName in ListRNAs:
        for target in Experiments:
            fileName = sampleName+"."+target+"_"+sampleName+".txt"

            if fileName in folderFileList:
                filePath = os.path.join(mutationstringDir,fileName)
                txtFileList.append(fileName)
                refPath =  os.path.join(fastaDir,sampleName+"."+FastaFileExtension)
                csvDir=os.path.join(outputDir, "counted_mutations")
                Fct.CheckFolderExist(csvDir)
                outPath = os.path.join(csvDir,sampleName+"_"+target+".csv")
                # mutations with associated phred scores below the minPhredToCount threshold are ignored
                argList = ["./countMutations",
                            "-file_in", filePath,
                            "-ref_seqs", refPath,
                            "-sample_name", sampleName,
                            "-target_name", sampleName,
                            "-file_out", outPath,
                            "-min_phred", str(minphredtocount)]

                argListDict[fileName] = argList


    stdOuts, stdErrs = Fct.spawnProcesses(maxProc, txtFileList, argListDict)

    for fileName in stdErrs:

        stdErr = stdErrs[fileName]
        stdOut = stdOuts[fileName]
        if stdOut != "" and stdOut is not None:
            print str(stdOut)+"\n"
        if len(stdErr) > 0:
            errorString = "Error: file %s mutation counting failed.\n"%(fileName)
            errorString += "\n".join(stdErr.splitlines())
            print errorString+"\n"

            print errorString
            sys.exit(1)
        else:
            print "File %s mutations counted successfully.\n"%(fileName)

            print 'mutations counted successfully'
except:
    errorString = "Error:" + traceback.format_exc() + "\n"
    errorString += "Mutation counting failed %s."% Fct.timeStamp()
    print errorString+"\n"

    print errorString
    sys.exit(1)
'''
#---------------------------------------------------------------------------------------
#3 Generate final reactivity profiles

print ("\nStarting reactivity profile creation at %s\n"%(Fct.timeStamp()))


#csvDir = os.path.join(outputDir,"counted_mutations/")
Fct.CheckFolderExist(os.path.join(outputDir, "reactivity_profiles"))
try:
    argListDict = {}
    csvDir = os.path.join(outputDir, "counted_mutations/")
    for profileName in ListRNAs:
        argList = ["python",
                    './generateReactivityProfiles.py',
                    csvDir,
                    profileName,
                    profileName,
		            csvDir+profileName+"_"+config.get("Parameters","Experimentshape")+".csv",
                    csvDir+profileName+"_"+config.get("Parameters","Experimentuntreated")+".csv",
            	    csvDir+profileName+"_"+config.get("Parameters","Experimentdenatured")+".csv"]

        argListDict[profileName] = argList



    stdOuts, stdErrs = Fct.spawnProcesses( maxProc,ListRNAs, argListDict)

    for profileName in stdErrs:
        stdErr = stdErrs[profileName]
        stdOut = stdOuts[profileName]
        if stdOut != "" and stdOut is not None:
            print stdOut+"\n"
        if "Error" in stdErr:
            errorString = "Error: reactivity profile %s generation failed.\n"%(profileName)
            errorString += "\n".join(stdErr.splitlines())
            print (errorString+"\n")

            print errorString
            sys.exit(1)
        else:
            print("Reactivity profile %s generated successfully.\n"%(profileName))


except:
    errorString = "Error:" + traceback.format_exc() + "\n"
    errorString += "Reactivity profile creation failed %s."%Fct.timeStamp()
    print(errorString+"\n")

    print errorString
    sys.exit(1)

doneMessage = "Pipeline completed %s."%Fct.timeStamp()
print doneMessage+"\n"
