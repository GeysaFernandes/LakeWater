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

import numpy as np
from numpy import trapz

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from sklearn.decomposition import PCA
from sklearn.neighbors import KDTree

import generate_fake_lake as fl

%matplotlib inline


# ****** functions ********

#read n .fna database files in the specified path
#set n = 0 to read all files
def ReadDataBaseFilenames(_path, n, filename_filename):
    seqList = []
    filenameList = []
    from os import path
    files = os.listdir(_path) #makes a list of all files in folder
    i = 0
    j = 0
    if '.DS_Store' in files:
        files.remove('.DS_Store')
    for f in files:
        for seq_record in SeqIO.parse(_path + f, "fasta"):
            s = seq_record.seq
            seqList.append(s) # reads each file into a list
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
            arr[dictionary[w]] += 1
        except:
            i = i
        i += 1
        if(i+n > sLen):
            break

    return arr


def kdtree(data, lake_matrix, k_neighbors = 10, leaf_size = 20):
    # training
    # kdtree = KDTree(data, leaf_size=leaf_size, metric='minkowski')
    kdtree = KDTree(data, leaf_size=leaf_size, metric='euclidean')

    # testing
    distances, indices = kdtree.query(lake_matrix, k=k_neighbors)
    return np.array(indices), distances

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

def process():
    lake_matrix = []

    n = 4; D = CreateDictionary(n)
    for w in lake:
        arr = FeatureVector(D, str(w), n)
        arr = np.divide(np.array(arr), len(w))
        lake_matrix.append(arr)

    len_lake = len(lake)
    len_viruses = len(virus_matrix)

    matrix = np.vstack((virus_matrix,bact_matrix))
    virus_bact_filenames = np.hstack((virus_filenames,bact_filenames))

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
    siz1, siz2 = data.shape
    neighbors = round(perc*siz1)
    print("running knn with ", neighbors, " neighbors....")

    # classification call
    indices, dist = kdtree(data, lake_matrix, k_neighbors = neighbors, leaf_size = 30)

    return indices, dist

# ******************************************************* main ****************************************************

# to use preprocessed bacteria and virus files
use_presaved = True
presaved_path = "../presaved/"

# if not using presaved, enter here where are the .fasta files
virus_database_path = "../database/virus/"
bact_database_path = "../database/bact/"

lake_path = "../database/lake/"
lake_quantity = 1000

# to reduce dimensionality of the data
use_pca = False
pca_components = 100

# percentage of neighbors to look for
perc = .15

if use_presaved:
    virus_matrix = np.loadtxt(presaved_path + "virus_features.txt")
    bact_matrix = np.loadtxt(presaved_path + "bact_features.txt")
    virus_filenames = filenameToMatrix(presaved_path + "virus_filenames.txt")
    bact_filenames = filenameToMatrix(presaved_path + "bact_filenames.txt")
else:
    print("reading data...")
    known_viruses, virus_filenames = ReadDataBaseFilenames(virus_database_path, 0, presaved_path + "virus_filenames.txt")
    known_bacterias, bact_filenames = ReadDataBaseFilenames(bact_database_path, 0, presaved_path + "bact_filenames.txt")
    virus_matrix = []
    bact_matrix = []
    n = 4; D = CreateDictionary(n)
    print("generating viruses feature vectors...")
    for w in known_viruses:
        arr = FeatureVector(D, str(w), n)
        arr = np.divide(np.array(arr), len(w))
        virus_matrix.append(arr)
    m=0
    print("generating bacterias feature vectors...")
    for w in known_bacterias:
        arr = FeatureVector(D, str(w), n)
        arr = np.divide(np.array(arr), len(w))
        m+=1
        bact_matrix.append(arr)
    # saving feature vectors in file
    fid1 = open(presaved_path + 'virus_features.txt', 'w')
    fid2 = open(presaved_path + 'bact_features.txt', 'w')
    fid1.write(listToString(virus_matrix))
    fid2.write(listToString(bact_matrix))
    fid1.close(); fid2.close()


print("reading lake...")
lake, lake_filenames = ReadDataBaseFilenames(lake_path, lake_quantity, presaved_path + "lake_filenames.txt")
print(np.array(lake).shape)

indices, dist = process()

print("finished")
