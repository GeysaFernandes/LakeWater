#Alana Laryssa, Cauanne Linhares, Geysa Fernandes
#Brazilian Scientific Mobility Program 2015-2016

#Libraries import may not work dependind on the enviroment you are using, we are using the Jupyter Notebook
from random import choice
from string import ascii_uppercase
import math
import time
import numpy as np
from scipy.integrate import simps
from numpy import trapz
import random
import operator
import itertools
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#####################################################
# Creates a txt file called lake.txt with real data
# from bacterial and viral database
#
#####################################################


# ****** functions ********

#read n .fna database files in the specified path
#set n = 0 to read all files
def ReadDataBase(_path, n):
    seqList = []
    from os import path
    files = os.listdir(_path) #makes a list of all files in folder
    i = 0
    j = 0
    for f in files:
        for seq_record in SeqIO.parse(_path + f, "fasta"):
            seqList.append(seq_record.seq) # reads each file into a list
            j += 1
            if(n > 0):
                if(j > n-1):
                    i = n
                    break
        i += 1
        j = 0

        if(n > 0):
            if(i > n-1):
                break

    return seqList

#creates dictionary with all permutations of length n
#with repetition to index Feature Vector
def CreateDictionary(n):
    chars = "ACGT"
    arr = list(itertools.product(chars, repeat=n))

    D = {}
    i = 0

    for a in arr:
        D[''.join(a)] = i
        i += 1

    return D

#builds the feature vector for sequence using specified indexing dictionary
def FeatureVector(dictionary, sequence, n):
    sLen = len(sequence)
    arr = [0]*4**n
    i = 0

    while(1):
        w = sequence[i:i+n]
        try:
            arr[D[w]] += 1
        except:
            i = i
        i += 1
        if(i+n > sLen):
            break

    return arr

#Reads the DB files and puts the information of the file in a array of strings
def readfile(filename):
    temp = open(filename, 'r').read().split('\n')
    return temp


#returns a random string of specified length
#length: strign length
def randomword(length):
    return (''.join(choice('ACGT') for i in range(0, length)))

#retuns an array of random strings
#size: how many strings there will be in the array
#lakeMinLen: min sequence length
#lakeMaxLen: max sequence length
def lakeString(size, lakeMinLen, lakeMaxLen):
    lake_water = []
    for i in range(0, size):
        random.seed()
        #generates a random sequence length
        y = random.randint(lakeMinLen, lakeMaxLen)

        _str = randomword(y)
        lake_water.append(_str)
    return lake_water

#configures the controlled data
#parameters:
#k_v: known_virus
#k_b: known_bacterias
#nControlled: how many records will be controlled
#error: corruption limit
def dataController(k_v, k_b, nControlVir, nControlBact, lakeMinLen, lakeMaxLen, indexes):
    controlled_lake = []
    vLen = len(k_v)
    bLen = len(k_b)

    #creates the specified amount of controlled lake sequences with viral pieces
    for j in range (0, nControlVir):
        random.seed()
        #get a random virus
        ind = random.randint(1, vLen-1)
        indexes.append(ind)
        vir = k_v[ind]
        #get the length of the random virus
        virLen = len(vir)
        count = 0
        #tries to find valid random substrings to be used
        num1 = random.randint(10, virLen)
        num2 = num1 - random.randint(1, num1)   
        print("N1:", num1, "N2:", num2)    

        #get a substring
        strg = vir[num2:num1]
        controlled_lake.append(strg)

    #creates the specified amount of controlled lake sequences with bacterial pieces
    for j in range (0, nControlBact):
        random.seed()
        ind = random.randint(1, bLen-1)
        indexes.append(ind + vLen)
        bact = k_b[ind]
        bactLen = len(bact)
        count = 0
        #tries to find valid random substrings to be used
        num1 = random.randint(10, bactLen)
        num2 = num1 - random.randint(1, num1)    
        print("N1:", num1, "N2:", num2)   

        #get a substring
        strg = bact[num2:num1]
        controlled_lake.append(strg)        

    return controlled_lake


#configures the controlled data
#parameters:
#k_v: known_virus
#k_b: known_bacterias
#nControlled: how many records will be controlled
#error: corruption limit
def dataController2(k_v, k_b, nControlVir, nControlBact, lakeMinLen, lakeMaxLen, indexes):
    controlled_lake = []
    vLen = len(k_v)
    bLen = len(k_b)

    #creates the specified amount of controlled lake sequences with viral pieces
    for j in range (0, nControlVir):
        random.seed()
        print(j)
        #get a random virus
        ind = random.randint(1, vLen-1)
        indexes.append(k_v[ind])
        for seq_record in SeqIO.parse(k_v[ind], "fasta"):
            vir = seq_record.seq
        #get the length of the random virus
        virLen = len(vir)
        count = 0
        #tries to find valid random substrings to be used
        num1 = random.randint(10, virLen)
        num2 = num1 - random.randint(1, num1)    
        print("N1:", num1, "N2:", num2)   

        #get a substring
        strg = vir[num2:num1]
        controlled_lake.append(strg)

    #creates the specified amount of controlled lake sequences with bacterial pieces
    for j in range (0, nControlBact):
        random.seed()
        print(j)
        ind = random.randint(1, bLen-1)
        indexes.append(k_b[ind])
        for seq_record in SeqIO.parse(k_v[ind], "fasta"):
            bact = seq_record.seq
        bactLen = len(bact)
        count = 0
        #tries to find valid random substrings to be used
        num1 = random.randint(10, bactLen)
        num2 = num1 - random.randint(1, num1)    
        print("N1:", num1, "N2:", num2)

        #get a substring
        strg = bact[num2:num1]
        controlled_lake.append(strg)
        
    return controlled_lake, indexes

def Generate_lake(known_viruses, known_bacterias, nControlVir, nControlBact):
    indexes = []
    lakeMinLen = 500
    lakeMaxLen = 2000
    # lake = dataController(known_viruses, known_bacterias, nControlVir, nControlBact, lakeMinLen, lakeMaxLen, indexes)
    lake, lake_index = dataController2(known_viruses, known_bacterias, nControlVir, nControlBact, lakeMinLen, lakeMaxLen, indexes)

    #save in file
    fd = open("../database/lake.txt", "w+")
    llen = len(lake)
    fd.write(str(nControlVir)+'\n')
    fd.write(str(nControlBact)+'\n')
    for i in range(0, llen):
        if(i > 0):
            fd.write('\n')
        fd.write(str(lake[i]))
    fd.close()

    fid = open("../database/lake_filenames.txt", 'w+')
    for i in indexes:
        fid.write(i + '\n')
    fid.close()

    return indexes

# ******************************************************* main ****************************************************

#test = readfile("../database/lake.txt")
#print(test)
