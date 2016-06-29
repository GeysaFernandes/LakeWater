#Alana Laryssa, Cauanne Linhares, Geysa Fernandes
#Brazilian Scientific Mobility Program 2015-2016

#Libraries import may not work dependind on the enviroment you are using, we are using the Jupyter Notebook
from random import choice
from string import ascii_uppercase
import math
import time
from swalign import swalign #The swalign file must be downloaded and kept in the same folder as the .py codes
from scipy.stats import beta
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
from numpy import trapz
import random
import operator
import itertools
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import KMeans

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
        #get a random virus
        ind = random.randint(1, vLen-1)
        indexes.append(ind)
        vir = k_v[ind]
        #get the length of the random virus
        virLen = len(vir)            
        count = 0
        #tries to find valid random substrings to be used
        while(True):
            num1 = random.randint(0, virLen)
            num2 = num1 + random.randint(lakeMinLen, lakeMaxLen)
            if num2 > virLen-1:
                continue
            else:
                break
        if num1 < num2:
            #controls the max length
            if(num2 - num1) > lakeMaxLen:
                num2 = num1 + lakeMaxLen
                
            #get a substring
            strg = vir[num1:num2]
            controlled_lake.append(strg)
        elif num2 < num1:
            #controls the max length
            if(num1 - num2) > lakeMaxLen:
                num1 = num2 + lakeMaxLen
                
            #get a substring
            strg = vir[num2:num1]
            controlled_lake.append(strg)
        else:
            #controls the max length
            if(virLen - num2) > lakeMaxLen:
                num2 += lakeMaxLen
            else:
                num2 = virLen-1
                
            #get a substring
            strg = vir[num1:num2]
            controlled_lake.append(strg)  
    
    #creates the specified amount of controlled lake sequences with bacterial pieces
    for j in range (0, nControlBact):
        ind = random.randint(1, bLen-1)
        indexes.append(ind)
        bact = k_b[ind]
        bactLen = len(bact)            
        count = 0
        while(True):
            num1 = random.randint(0, bactLen)
            num2 = num1 + random.randint(lakeMinLen, lakeMaxLen)
            if num2 > bactLen-1:
                continue
            else:
                break
        if num1 < num2:
            #controls the max length
            if(num2 - num1) > lakeMaxLen:
                num2 = num1 + lakeMaxLen
                
            #get a substring
            strg = bact[num1:num2]
            controlled_lake.append(strg)
        elif num2 < num1:
            #controls the max length
            if(num1 - num2) > lakeMaxLen:
                num1 = num2 + lakeMaxLen
                
            #get a substring
            strg = bact[num2:num1]
            controlled_lake.append(strg)
        else:
            #controls the max length
            if(virLen - num2) > lakeMaxLen:
                num2 += lakeMaxLen
            else:
                num2 = bactLen-1
                
            #get a substring
            strg = bact[num1:num2]
            controlled_lake.append(strg)
    
    return controlled_lake

def Generate_lake(known_viruses, known_bacterias):
    print("generate lake")
    print(known_viruses)
    indexes = []
    print("before:", indexes)
    nControlVir = 20
    nControlBact = 20
    lakeMinLen = 200
    lakeMaxLen = 500
    lake = dataController(known_viruses, known_bacterias, nControlVir, nControlBact, lakeMinLen, lakeMaxLen, indexes)    
    print("after:", indexes)
    print("finished reading data...")
    print("writing file...")

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

    print("end of execution...")
    return indexes

def Teste(string):
    print(string)
# ******************************************************* main ****************************************************

#test = readfile("../database/lake.txt")
#print(test)