
# coding: utf-8

# In[3]:

from random import choice
from string import ascii_uppercase
import math
import time
from swalign import swalign
from scipy.stats import beta
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
from numpy import trapz
import random
import operator

print(swalign)
# ****** functions ********
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

#returns one array of strings that represents database data set
#dbMinSeqLen: min string length
#dbMaxSeqLen: max string length
#dbSize: max number sequences in the array
def dbStrings(dbMinSeqLen, dbMaxSeqLen, dbSize):
    _array = []
    #generates how many different sequences there will be for the array
    random.seed()
    #x = random.randint(2, dbSize)
    x = 30

    
    for i in range(0, x):
        #generates a random sequence size
        random.seed()
        y = random.randint(dbMinSeqLen, dbMaxSeqLen)
        _str = randomword(y)
        _array.append(_str)
    return _array

#corrupts the string if it is under the corruption rate
def stringCorruption(strg, corrupt_rate):
    _str = ''
    for k in strg:
        #get a random number between 0 and 1
        random.seed()
        floatr = random.random()
        if  floatr > corrupt_rate :
            _str = _str + k
        else:
            _str = _str + choice('ACGT')     
    return _str

#configures the controlled data
#parameters:
#k_v: known_virus
#k_b: known_bacterias
#nControlled: how many records will be controlled
#error: corruption limit
def dataController(k_v, k_b, error, nControlVir, nControlBact, lakeMinLen, lakeMaxLen):
    controlled_lake = []
    vLen = len(known_viruses)
    bLen = len(known_bacterias)
    
    #creates the specified amount of controlled lake sequences with viral pieces
    for j in range (0, nControlVir):
        #get a random virus
        vir = known_viruses[random.randint(1, vLen-1)]
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
            #corrupts the string if it's under the error rate
            strg = stringCorruption(strg, error)
            controlled_lake.append(strg)
        elif num2 < num1:
            #controls the max length
            if(num1 - num2) > lakeMaxLen:
                num1 = num2 + lakeMaxLen
                
            #get a substring
            strg = vir[num2:num1]
            #corrupts the string if it's under the error rate
            strg = stringCorruption(strg, error)
            controlled_lake.append(strg)
        else:
            #controls the max length
            if(virLen - num2) > lakeMaxLen:
                num2 += lakeMaxLen
            else:
                num2 = virLen-1
                
            #get a substring
            strg = vir[num1:num2]
            #corrupts the string if it's under the error rate
            strg = stringCorruption(strg, error)
            controlled_lake.append(strg)  
    
    #creates the specified amount of controlled lake sequences with bacterial pieces
    for j in range (0, nControlBact):
        bact = known_bacterias[random.randint(1, bLen-1)]
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
            #corrupts the string if it's under the error rate
            strg = stringCorruption(strg, error)
            controlled_lake.append(strg)
        elif num2 < num1:
            #controls the max length
            if(num1 - num2) > lakeMaxLen:
                num1 = num2 + lakeMaxLen
                
            #get a substring
            strg = bact[num2:num1]
            #corrupts the string if it's under the error rate
            strg = stringCorruption(strg, error)
            controlled_lake.append(strg)
        else:
            #controls the max length
            if(virLen - num2) > lakeMaxLen:
                num2 += lakeMaxLen
            else:
                num2 = bactLen-1
                
            #get a substring
            strg = bact[num1:num2]
            #corrupts the string if it's under the error rate
            strg = stringCorruption(strg, error)
            controlled_lake.append(strg)
    
    return controlled_lake
    
# ******************************************************* main ****************************************************

####---------------- config. variables --------------------------------####
percentage = 0      #how many will have virus/bacteria pieces #0 for no control
bacterialRate = 0.3 #how many of the controlled lake sequences will bacterial piece
error = 0.0         #corruption rate
distLimit = 0       #distance between matching sequences

lakeSize = 10    #max number of lake sequences
lakeMinLen = 10     #min lake string length
lakeMaxLen = 100     #max lake string length

dbSize = 15  #max number sequences in the array
dbMaxSeqLen = 200  #max db string length
dbMinSeqLen = lakeMinLen #min db string length

#print configuration
print("\n\nConfig. Variables")
print("  Controlled Data:", percentage*100, "%")
print("  Percentage of Bact. data:", bacterialRate*100, "%")
print("  Corruption Rate:", error*100, "%")
print("  Ditance Limit:", distLimit)
print("  Db string min lenght:", dbMinSeqLen)
print("  Db string max lenght:", dbMaxSeqLen)
print("  Db max size:", dbSize)
print("  Lake string min lenght:", lakeMinLen)
print("  Lake string max lenght:", lakeMaxLen)
print("  Lake data max size:", lakeSize, "\n\n")

#generates how many different sequences there will be for the lake
#random.seed()
#size = random.randint(2, lakeSize)
size = 10

#generates the database representation
known_viruses = dbStrings(dbMinSeqLen, dbMaxSeqLen, dbSize)
known_bacterias = dbStrings(dbMinSeqLen, dbMaxSeqLen, dbSize)

#generates the lake sequences w/ the specified percentage of controlled data
P = math.ceil(size * percentage)
lake = lakeString(size - P, lakeMinLen, lakeMaxLen) #produce a part with totally random sequences  

#calculates how many sequences will be viral and how many will be bacterial
nControlBact = int(math.ceil(P * bacterialRate))
if bacterialRate == 0.0:
    nControlVir = P
else:
    nControlVir = int((nControlBact * (1 - bacterialRate)) / bacterialRate)

print("# inserted viral pieces:", nControlVir, "\n# inserted bacterial pieces: ", nControlBact)
size = size - P + nControlVir + nControlBact
print("True bacterial percentage:", nControlBact/size*100, "%")

#produce the other part with viral or bacterial pieces
controlled = dataController(known_viruses, known_bacterias, error, nControlVir, nControlBact, lakeMinLen, lakeMaxLen)
lake = lake + controlled #join the sequences

#Lakewater String

#For 100% controlled
fd = open("lakewater_string_not_controlled.txt", "w+")
for i in lake:
    fd.write(i)
    fd.write('\n')
fd.close()

#Virus String
fd_v = open("virus_string_not_controlled.txt", "w+")
for i in known_viruses:
    fd_v.write(i)
    fd_v.write("\n")
fd_v.close()

#Bacteria String
fd_b = open("bacteria_string_not_controlled.txt", "w+")
for i in known_bacterias:
    fd_b.write(i)
    fd_b.write("\n")
fd_b.close()


# In[ ]:




# In[ ]:



