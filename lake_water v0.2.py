
# coding: utf-8

# In[13]:



# In[16]:

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

#returns one array of strings that represents database data set
#dbMinSeqLen: min string length
#dbMaxSeqLen: max string length
#dbSize: max number sequences in the array
def dbStrings(dbMinSeqLen, dbMaxSeqLen, dbSize):
    _array = []
    #generates how many different sequences there will be for the array
    random.seed()
    x = random.randint(2, dbSize)

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
percentage = 1      #how many will have virus/bacteria pieces
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
random.seed()
size = random.randint(2, lakeSize)

#generates the database representation
#known_viruses = dbStrings(dbMinSeqLen, dbMaxSeqLen, dbSize)
#known_bacterias = dbStrings(dbMinSeqLen, dbMaxSeqLen, dbSize)
known_viruses = readfile("virus_string_controlled.txt")
known_bacterias = readfile("bacteria_string_controlled.txt")

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
#lake = lake + controlled #join the sequences
lake = readfile("lakewater_string_conntrolled.txt")

#testing sequences
#print(lake)
#print("\n")

#print(known_bacterias)
#print("\n")

#print(known_viruses)
#print("\n")

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

start = time.clock() #start clock to measure proccessing time

for l in lake:       #for each sequece in the lake
    #print("***l len: ", len(l))
    viralMatches = 0 
    bactMatches = 0
    lLen = len(l)                                  #get sequence length
    #gives different weights depending of the sequence length 
    pnts = math.ceil((lLen/dbMaxSeqLen*100) + 1)
    _limit = round(lLen * (100 - distLimit) / 100) #get the minimum points to be considered a match
    for v in known_viruses:                        #compare the sequence with all virus in the database
        alignment = sw.align(v, l)
        if alignment.score >= _limit:
            viralMatches += pnts       
    for b in known_bacterias: #compare the sequence with all bacterias in the database
        alignment = sw.align(b, l)
        if alignment.score >= _limit:
            bactMatches += pnts

    #draws beta distributions    
    bactMatches += 1
    viralMatches += 1    
    mean, var, skew, kurt = beta.stats(bactMatches, viralMatches, moments='mvsk')
    #_label = "Viral:", viralMatches, "Bacterial:", bactMatches
    #ax.plot(x, beta.pdf(x, bactMatches, viralMatches),'p-', lw=3, alpha=0.6, label=_label)
    y = beta.pdf(x, bactMatches, viralMatches)
    yIndex = np.argmax(_y)

    #multiplies the coordinates
    _y *= y 

    
end = time.clock()
print("\n\n>>>Elapsed Time:", round(end - start, 3))    

#Gets the highest value, the peak value
yIndex = np.argmax(_y)

print("\n\nReport")
print("   Viruses:", round((1 - x[yIndex]) * 100, 2), "%")
print("   Bacteria:", round(x[yIndex] * 100, 2), "%")

#plots the multiplication between the beta distribution curves 
fig, ax2 = plt.subplots(1, 1)
ax2.plot(x, _y, 'y-', lw = 5, alpha=0.9, label="label")
ax2.set_title("Beta Distribution Multiplication")
ax2.set_xlabel("Bacterial Rate")
ax2.legend(loc='best', frameon=False)
plt.show()



# In[ ]:



