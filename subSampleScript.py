#!/usr/bin/env python
import random

def subSample(listOfFileNames):
    for fileName in listOfFileNames:
        file = open(fileName, 'r')
        lines = file.readlines()
        randomLines = random.sample(lines, 10)
        file.close()

        newFile = open(fileName + '.subSample', 'w')
        [ newFile.write(l) for l in randomLines ]
        newFile.close()

listOfFileNames = ["run_files/RAxML_bootstrap.orfg1.last_2", "run_files/RAxML_bootstrap.orfg10.last_2", "run_files/RAxML_bootstrap.orfg10_5.last_3", "run_files/RAxML_bootstrap.orfg11.last_2", "run_files/RAxML_bootstrap.orfg12.last_2", "run_files/RAxML_bootstrap.orfg13.last_2", "run_files/RAxML_bootstrap.orfg14.last_2", "run_files/RAxML_bootstrap.orfg15.last_2", "run_files/RAxML_bootstrap.orfg2.last_2", "run_files/RAxML_bootstrap.orfg3.last_2", "run_files/RAxML_bootstrap.orfg3_5.last_2", "run_files/RAxML_bootstrap.orfg4.last_2", "run_files/RAxML_bootstrap.orfg5.last_2", "run_files/RAxML_bootstrap.orfg6.last_2", "run_files/RAxML_bootstrap.orfg7.last_2", "run_files/RAxML_bootstrap.orfg8.last_2", "run_files/RAxML_bootstrap.orfg9.last_2"]

subSample(listOfFileNames)
