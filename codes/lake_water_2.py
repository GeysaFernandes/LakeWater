# Adding this snippet so the code can run on osx
from sys import platform as platform_name
if platform_name == "darwin":
   import sys
   sys.path.append('//anaconda/lib/python3.5/site-packages/')

import random
from random import choice

import operator
import os
import math
import time
from string import ascii_uppercase

from swalign import swalign

import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from numpy import trapz

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import KDTree

import threading

import generate_fake_lake as fl


# ****** aux classes ********
class featureVectorThread (threading.Thread):
    def __init__(self, matrix, i, D, w, n):
        threading.Thread.__init__(self)
        self.threadID = i
        self.dict = D
        self.sequence = w
        self.n = n
        self.matrix = matrix
    def run(self):
        print("Starting ", self.i)
        arr = FeatureVector(D, str(w), n)
        arr = np.divide(np.array(arr), len(w))
        matrix[i] = arr
        print("Exiting ", self.i)

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


#read n .fna database files in the specified path
#set n = 0 to read all files
def ReadDataBaseFilenames(_path, n, filename_filename):
    seqList = []
    filenameList = []
    from os import path
    files = os.listdir(_path) #makes a list of all files in folder
    i = 0
    j = 0
    for f in files:
        for seq_record in SeqIO.parse(_path + f, "fasta"):
            seqList.append(seq_record.seq) # reads each file into a list
            filenameList.append(_path + f)
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

    fid = open(filename_filename, 'w')
    for f in filenameList:
        fid.write(f + '\n')
    fid.close()
    return seqList, filenameList


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


def kdtree(data, lake_matrix, k_neighbors = 10, leaf_size = 20):
    # training
    # kdtree = KDTree(data, leaf_size=leaf_size, metric='minkowski')
    kdtree = KDTree(data, leaf_size=leaf_size, metric='euclidean')

    # testing
    distances, indices = kdtree.query(lake_matrix, k=k_neighbors)
    return np.array(indices), distances


def clustering(data, lake_matrix, num_clusters = 12):
    # training
    estimator = KMeans(n_clusters=num_clusters)
    estimator.fit(data)
    training_labels = estimator.labels_

    clusters = [[]]*num_clusters  # stores the index of the points in each cluster
    for clust in range(0,num_clusters):
        labels_idx = np.where(training_labels == clust)[0]
        clusters[clust] = labels_idx

#     # printing cluster data
#     i = 0
#     for clt in clusters:
#         virs = [c for c in clt if c < len_viruses]
#         perc = len(virs)/len(clt)
#         print("% of virus in cluster ", i, ": ", perc)
#         i += 1

    # testing
    indices = []
    estimated_labels = estimator.predict(lake_matrix)
    for lbl in estimated_labels:
        indices.append(clusters[lbl])

    return np.array(indices)

def check_indices(indices, lake_indices, filenames):
    right = 0
    cnt = 0

    rightv = 0
    rightb = 0

    for i in lake_indices:
        # print(i, " -- ", indices[cnt])
        for ii in indices[cnt]:
            if i == filenames[ii]:
                if cnt < 10:
                    rightv += 1
                else:
                    rightb += 1
                # print("right ", cnt)
                right += 1
                # print("achou em: ", ii)
                break
        cnt += 1
    print ("right answers: ", right, " / ", len(lake_indices))
    print ("bacteria: ", rightb, " | virus: ", rightv)
    return (right / len(lake_indices))

def listToString(list):
    string = ''
    for ll in list:
        for l in ll:
            string += str(l) + ' '
        string += '\n'
    return string

def filenameToMatrix(file):
    l = []
    fid = open(file,'r')
    for f in fid:
        l.append(f.replace('\n',''))
    return l

def computeFeatureVector(known_bacterias):
    import multiprocessing as mp
    num_thread = multiprocessing.cpu_count()
    matrix = [[]]*len(known_bacterias)
    for i in range(num_thread):
        featureVectorThread(matrix, i, D, w, n)

    n = 4; D = CreateDictionary(n)
    for w in known_bacterias:
        arr = FeatureVector(D, str(w), n)
        arr = np.divide(np.array(arr), len(w))
        m+=1
        if(m%50 == 0):
            print(m,"of", len(known_bacterias))
        bact_matrix.append(arr)


# ******************************************************* main ****************************************************

use_presaved = True

if use_presaved:
    virus_matrix = np.loadtxt('../presaved/virus_features.txt')
    bact_matrix = np.loadtxt('../presaved/bact_features.txt')
    virus_filenames = filenameToMatrix("../presaved/virus_filenames.txt")
    bact_filenames = filenameToMatrix("../presaved/bact_filenames.txt")
else:
    print("reading data...")
    known_viruses, virus_filenames = ReadDataBaseFilenames("../database/virus/", 0, "../presaved/virus_filenames.txt")
    known_bacterias, bact_filenames = ReadDataBaseFilenames("../database/bact/", 0, "../presaved/bact_filenames.txt")
    virus_matrix = []
    bact_matrix = []
    n = 4; D = CreateDictionary(n)
    print("generating feature vector virus...")
    for w in known_viruses:
        arr = FeatureVector(D, str(w), n)
        arr = np.divide(np.array(arr), len(w))
        virus_matrix.append(arr)
    m=0
    print("generating feature vector bacteria...")
    for w in known_bacterias:
        arr = FeatureVector(D, str(w), n)
        arr = np.divide(np.array(arr), len(w))
        m+=1
        if(m%50 == 0):
            print(m,"of", len(known_bacterias))
        bact_matrix.append(arr)
    # saving feature vectors in file
    fid1 = open('../presaved/virus_features.txt', 'w')
    fid2 = open('../presaved/bact_features.txt', 'w')
    fid1.write(listToString(virus_matrix))
    fid2.write(listToString(bact_matrix))
    fid1.close(); fid2.close()


use_pca = True
pca_components = 50

x_data = []
y_data = []

# for it in range(1,25,1):
for it in [1,3,5,10]:

    x_data.append(it)
    print("===================== ", it)

    # lake = ReadDataBase("../database/lake/", 50)

    print("generating lake...")
    lake_indices = fl.Generate_lake(virus_filenames, bact_filenames, 10, 10)
    lake = readfile("../database/lake.txt")
    num_virus_lake = int(lake[0])
    num_bact_lake = int(lake[1])
    lake = lake[2:]

    print("finished reading data")

    lake_matrix = []

    n = 4; D = CreateDictionary(n)
    for w in lake:
        arr = FeatureVector(D, str(w), n)
        arr = np.divide(np.array(arr), len(w))
        lake_matrix.append(arr)

    len_lake = len(lake)
    len_viruses = len(virus_matrix)

    print("lake shape: ", np.matrix(lake_matrix).shape)
    print("virus shape: ", np.matrix(virus_matrix).shape)
    print("bacteria shape: ", np.matrix(bact_matrix).shape)

    matrix = np.vstack((virus_matrix,bact_matrix))

    ###### PCA
    if use_pca:
        print("performing PCA...")
        X = np.array(matrix)
        # PCA input: samples x features
        pca = PCA(n_components=pca_components)
        Xhat = pca.fit_transform(X)
        print("Percentage of represented variance: ", sum(pca.explained_variance_ratio_))


    ###### CLASSIFICATION
    if use_pca:
        data = np.array(Xhat)
        lake_matrix = pca.transform(lake_matrix)
    else:
        data = np.array(matrix)


    ### choose a classification method
    ## the method will pick the best candidates to perform local alignment in each lake sample
    perc = it/100
    neighbors = round(perc*len(matrix))
    print("running knn with ", neighbors, " neighbors....")
    # neighbors = len(matrix)
    indices, dist = kdtree(data, lake_matrix, k_neighbors = neighbors, leaf_size = 30)
    # indices = clustering(data, lake_matrix, num_clusters = 12)

    virus_bact_filenames = np.hstack((virus_filenames,bact_filenames))

    # print("lake atual: ", lake_indices[6])
    # print("dist:")
    # for i in range(0,100):
    #     print(indices[6,i], " | ", dist[6,i])

    # for l in lake_indices:
    #     print(l)

    acc = check_indices(indices, lake_indices, virus_bact_filenames)
    y_data.append(acc)

print("x: ", x_data)
print("y: ", y_data)

plt.plot(x_data, y_data)
plt.savefig('graph')
plt.show()

# ###### LOCAL ALIGNMENT
# match = 1                                                  # scores 1 point for matching letters
# mismatch = -1                                              # looses 1 point for mismatching letters
# scoring = swalign.NucleotideScoringMatrix(match, mismatch)
# sw = swalign.LocalAlignment(scoring)                       # you can also choose gap penalties, etc...
# beta_x = np.arange(0, 1, 0.001)  # for beta distribution
# _y = 1  # final curve
# idx = 0
# right = 0

# print("performing alingment...")
# start = time.clock()
# for lake_sample in lake:
#     ## setup for local alignment
#     viral_matches = 0
#     bact_matches = 0
#     sample_length = len(lake_sample)

#     ## temporary databases
#     virus_temp_db = []
#     bact_temp_db = []
#     print("performing with ", len(indices[idx]), " samples")
#     for i in indices[idx]:
#         if i < len_viruses:
#             virus_temp_db.append(known_viruses[i])
#         else:
#             bact_temp_db.append(known_bacterias[i-len_viruses])

#     if not(virus_temp_db):  # if virus is empty
#         print("   all samples are bacteria: ", idx)
#         bact_matches = 1
#         viral_matches = 0
#         virus_temp_db = []
#         bact_temp_db = []
#     if not(bact_temp_db):  # if bacteria is empty
#         print("   all samples are viruses: ", idx)
#         bact_matches = 0
#         viral_matches = 1
#         virus_temp_db = []
#         bact_temp_db = []

#     ## local alignment
#     for v in virus_temp_db:
#         alignment = sw.align(v, lake_sample)
#         viral_matches += alignment.score
#     for b in bact_temp_db:
#         alignment = sw.align(b, lake_sample)
#         bact_matches += alignment.score

#     # checking answers
#     if viral_matches > bact_matches and idx < len_viruses:
#         right += 1
#     if bact_matches > viral_matches and idx >= len_viruses:
#         right += 1

#     print(viral_matches, " / ", bact_matches)

#     idx += 1

# end = time.clock()
# print("\n\n>>>Elapsed Time:", round(end - start, 3))
# print("accuracy: ", right/(len_lake))

# plt.savefig('myfig')
