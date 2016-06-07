
# coding: utf-8

# In[1]:

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


# ****** functions ********
#returns a random string of specified length
def randomword(length):
    return (''.join(choice('ACGT') for i in range(0, length)))

#retuns an array of random strings
def lakeString(size, maxSeqLen):     
    lake_water = []
    for i in range(0, size):
        random.seed()
        #generates a random sequence size
        y = random.randint(5, maxSeqLen)
        _str = randomword(y)
        lake_water.append(_str)
    return lake_water

#returns one array of strings that represents database data set
def dbStrings(dbMaxSeqLen):
    _array = []
    #generates how many different sequences there will be for the array
    random.seed()
    x = random.randint(2, 10)

    for i in range(0, x):
        #generates a random sequence size
        random.seed()
        y = random.randint(5, dbMaxSeqLen)
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
def dataController(k_v, k_b, nControlled, error, bacterialRate, maxSeqLen):
    controlled_lake = []
    vLen = len(known_viruses)
    bLen = len(known_bacterias)
    nControlBact = round(nControlled * bacterialRate)
    nControlVir = nControlled - nControlBact
    print("# inserted viral pieces:", nControlVir, "\n# inserted bacterial pieces: ", nControlBact)
    
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
            num2 = random.randint(0, virLen)
            #min string length is 4 
            if((num1 - num2) > 3 or (num2 - num1) > 3):
                break
            #loop for max 3 times or get the first 4 nucleic
            if(count > 2):
                num1 = 0
                num2 = 4
                break                    
            count += 1

        if num1 < num2:
            #controls the max length
            if(num2 - num1) > maxSeqLen:
                num2 = num1 + maxSeqLen
                
            #get a substring
            strg = vir[num1:num2]
            #corrupts the string if it's under the error rate
            strg = stringCorruption(strg, error)
            controlled_lake.append(strg)
        elif num2 < num1:
            #controls the max length
            if(num1 - num2) > maxSeqLen:
                num1 = num2 + maxSeqLen
                
            #get a substring
            strg = vir[num2:num1]
            #corrupts the string if it's under the error rate
            strg = stringCorruption(strg, error)
            controlled_lake.append(strg)
        else:
            #controls the max length
            if(virLen - num2) > maxSeqLen:
                num2 += maxSeqLen
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
            num2 = random.randint(0, bactLen)
            #min string length is 4 
            if((num1 - num2) > 3 or (num2 - num1) > 3):
                break
            # try 3 times to find valid random numbers or get the first 4 nucleic
            if(count > 2):
                num1 = 0
                num2 = 4
                break                    
            count += 1

        if num1 < num2:
            #controls the max length
            if(num2 - num1) > maxSeqLen:
                num2 = num1 + maxSeqLen
                
            #get a substring
            strg = bact[num1:num2]
            #corrupts the string if it's under the error rate
            strg = stringCorruption(strg, error)
            controlled_lake.append(strg)
        elif num2 < num1:
            #controls the max length
            if(num1 - num2) > maxSeqLen:
                num1 = num2 + maxSeqLen
                
            #get a substring
            strg = bact[num2:num1]
            #corrupts the string if it's under the error rate
            strg = stringCorruption(strg, error)
            controlled_lake.append(strg)
        else:
            #controls the max length
            if(virLen - num2) > maxSeqLen:
                num2 += maxSeqLen
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
percentage = 1      #how many will have virus/bacteria pieces
bacterialRate = 0.3 #from the controlled data how much will be bacterial data
error = 0.0         #corruption rate
distLimit = 0       #distance between matching sequences
maxLakeLen = 100    #max number of lake sequences
maxSeqLen = 150     #max lake string length
dbMaxSeqLen = 200  #max db string length

#print configuration
print("Config. Variables")
print("  Controlled Data:", percentage*100, "%")
print("  Corruption Rate:", error*100, "%")
print("  Ditance Limit:", distLimit)
print("  Lake string max lenght:", maxSeqLen, "\n\n")

#generates how many different sequences there will be for the lake
random.seed()
size = random.randint(2, maxLakeLen)

#generates the database representation
known_viruses = dbStrings(dbMaxSeqLen)
known_bacterias = dbStrings(dbMaxSeqLen)

#generates the lake sequences w/ the specified percentage of controlled data
P = math.ceil(size * percentage)
lake = lakeString(size - P, maxSeqLen)#produce a part with totally random sequences
controlled = dataController(known_viruses, known_bacterias, P, error, bacterialRate, maxSeqLen)#produce the other part with viral or bacterial pieces
lake = lake + controlled #join the sequences
i = len(lake)

#print sequences
#i = 0
#print("*** Lake Sequences")
#for w in lake:
#    print(i, " ", w)
#    i = i + 1 
#    
#print("\n\n*** Viral Sequences")
#m = 0
#for w in known_viruses:
#    print(m, " ", w)
#    m = m + 1 
#    
#print("\n\n*** Baterial Sequences")
#m = 0
#for w in known_bacterias:
#    print(m, " ", w)
#    m = m + 1 
#print("\n\n")

nBact = 0    #final counter
nVirus = 0   #final counter
nUnknown = 0 #final counter
_y = 1 #final curve

#Local Alignment settings
# choose your own values here… 2 and -1 are common.
match = 1                                                  #scores 1 point for matching letters
mismatch = -1                                              #looses 1 point for mismatching letters
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
sw = swalign.LocalAlignment(scoring)                       # you can also choose gap penalties, etc...

#draws beta distributions
x = np.arange(0, 1, 0.001)
fig, ax = plt.subplots(1, 1)

#start = time.clock() #start clock to measure proccessing time

for l in lake:       #for each sequece in the lake
    viralMatches = 0 
    bactMatches = 0
    lLen = len(l)                                  #get sequence length
    #gives different weights depending of the sequence length 
    pnts = round(lLen/dbMaxSeqLen*100)
    _limit = round(lLen * (100 - distLimit) / 100) #get the minimum points to be considered a match
    for v in known_viruses:                        #compare the sequence with all the virus in the database
        alignment = sw.align(v, l)
        if alignment.score >= _limit:
            viralMatches += pnts       
    for b in known_bacterias: #compare the sequence with all the virus in the database
        alignment = sw.align(b, l)
        if alignment.score >= _limit:
            bactMatches += pnts
    
    if bactMatches > viralMatches:
        nBact = nBact + 1
    elif viralMatches > bactMatches:
        nVirus = nVirus + 1 
    elif viralMatches == bactMatches and bactMatches > 0:
        nVirus = nVirus + 0.5
        nBact = nBact + 0.5
    else:
        nUnknown += 1        
    
    #draws beta distributions    
    bactMatches += 1
    viralMatches += 1    
    mean, var, skew, kurt = beta.stats(bactMatches, viralMatches, moments='mvsk')
    _label = "Viral:", viralMatches, "Bacterial:", bactMatches
    ax.plot(x, beta.pdf(x, bactMatches, viralMatches),'p-', lw=3, alpha=0.6, label=_label)
    y = beta.pdf(x, bactMatches, viralMatches)
    _y *= y    
        
#end = time.clock()
#print("\n\n>>>Elapsed Time:", round(end - start, 3))

print("\n\nReport")
print("   # Viral:", nVirus)
print("   # Bacterial:", nBact)
print("   # Unknown:", nUnknown)
print("   Viruses:", round(nVirus/i*100, 2), "%")
print("   Bacteria:", round(nBact/i*100, 2), "%")
print("   Unknown:", round(nUnknown/i*100, 2), "%")

#plots the multiplication between the beta distribution curves 
#shows all the figures  
fig, ax2 = plt.subplots(1, 1)
ax2.plot(x, _y, 'y-', lw = 5, alpha=0.9, label="label")
ax2.set_title("Beta Distribution Multiplication")
ax2.set_xlabel("Bacterial Rate")
#ax2.legend(loc='best', frameon=False) 

ax.set_title("Beta Distributions") 
ax.set_xlabel("Bacterial Rate")
#ax.legend(loc='best', frameon=False)  
plt.show()

