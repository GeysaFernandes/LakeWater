# Adding this snippet so the code can run on osx
from sys import platform as platform_name
if platform_name == "darwin":
   import sys
   sys.path.append('//anaconda/lib/python3.5/site-packages/')

from random import choice

import os
import math
import itertools
import time
from string import ascii_uppercase

import numpy as np

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from sklearn.decomposition import PCA
from sklearn.neighbors import KDTree

import generate_fake_lake as fl

# ****** functions ********

#read n .fna database files in the specified path
#set n = 0 to read all files
def ReadDataBaseFilenames(_path, n, filename_filename, species_filename, extension):
    seqList = []
    filenameList = []
    idList = []
    speciesList = []

    from os import path
    files = os.listdir(_path) #makes a list of all files in folder

    i = 0
    j = 0
    if '.DS_Store' in files:
        files.remove('.DS_Store')
    for f in files:
        for seq_record in SeqIO.parse(_path + f, extension):
            s = seq_record.seq
            seqList.append(s) # reads each file into a list
            filenameList.append(_path + f)
            idList.append(seq_record.id)

            desc = seq_record.description
            desc = desc.split('|')[-1]  # get the last entry
            desc = desc.split(',')[0]
            speciesList.append(desc)

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

    fid1 = open(species_filename, 'w')
    for f in speciesList:
        fid1.write(f + '\n')
    fid1.close()



    return seqList, filenameList, idList, speciesList


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

def getIDfromFilename(filename):
    return os.path.basename(filename).split('.')[0]

def process():
    lake_matrix = []

    nucleotides = 4; D = CreateDictionary(nucleotides)
    for w in lake:
        arr = FeatureVector(D, str(w), nucleotides)
        arr = np.divide(np.array(arr), len(w))
        lake_matrix.append(arr)

    len_lake = len(lake)
    len_viruses = len(virus_matrix)

    ###### PCA
    if use_pca:
        print("Reducing dimensionality from", matrix.shape[1], "to", pca_components, "...")
        X = np.array(matrix)
        # PCA input: samples x features
        pca = PCA(n_components=pca_components)
        Xhat = pca.fit_transform(X)
        print("Percentage of represented variance:", sum(pca.explained_variance_ratio_)*100)

    ###### CLASSIFICATION
    if use_pca:
        data = np.array(Xhat)
        lake_matrix = pca.transform(lake_matrix)
    else:
        data = np.array(matrix)

    ## the method will pick the best candidates to perform local alignment in each lake sample
    siz1, siz2 = data.shape
    neighbors = round(perc*siz1)
    print("Running k-NN with", neighbors, "neighbors....")

    # classification call
    indices, dist = kdtree(data, lake_matrix, k_neighbors = neighbors, leaf_size = 30)

    return indices, dist

def saveOutput():
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for i in range(0,len(lake)):
        curr = time.strftime("%d-%m-%Y--%H-%M-%S")
        fid = open(output_folder + curr + "--" + str(i), 'w+')

        for j in range(0,len(indices[i])):
            idx = indices[i][j]
            idx_str = getIDfromFilename(virus_bact_filenames[idx])
            spec = virus_bact_species[idx]
            line = ""
            line += ">[" + lake_ids[i] + "]"
            line += idx_str + "|" + spec + "|" + str(dist[i][j]) + "|"
            line += str(len(lake)) + "|" + str(len(indices[i]))
            fid.write(line + "\n")
            fid.write(str(lake[i]))
            fid.write("\n")
        fid.close()


# ******************************************************* parameters ****************************************************

# Set here if you want to use precomputed files.
# If yes, tell in what path the files are.
# In order to run the precomputed files,
# you have to run at least once without them.
use_presaved = True
presaved_path = "../presaved/"

# If you are not using precomputed files, enter below
# the name of the folders containing the .fasta files.
virus_database_path = "../database/virus/"
bact_database_path = "../database/bact/"
virus_extension = "fasta"
bact_extension = "fasta"

# Enter here the path to the folder with lake files.
# Below, enter how many lake samples you want to read.
lake_path = "../database/lake/"
lake_extension = "fastq"
lake_quantity = 5

# Set below the length of nucleotides sequence you want to
# consider when generating the feature vector
# Eg.: With 4 nucleotides, each position of the feature vector
# will count how many subsequences of 4 nucleotides are in the
# bacteria or virus sequence.
nucleotides = 4

# Set below if you want to reduce dimensionality of the data.
# Eg.: if you use 4 nucleotides, the feature vectores will be a
# 256-position vector. If you want to reduce the number of dimensions
# of this vector, insert below.
use_pca = False
pca_components = 100

# Enter the percentage of the total training data number that
# you want to choose as candidates
perc = .02

# Insert where you want to save the output files
output_folder = "../output/test0/"

# ******************************************************* main ****************************************************

if use_presaved:
    print("Reading presaved training data...")
    virus_matrix = np.loadtxt(presaved_path + "virus_features.txt")
    bact_matrix = np.loadtxt(presaved_path + "bact_features.txt")
    virus_filenames = filenameToMatrix(presaved_path + "virus_filenames.txt")
    bact_filenames = filenameToMatrix(presaved_path + "bact_filenames.txt")
    virus_species = filenameToMatrix(presaved_path + "virus_filenames.txt")
    bact_species = filenameToMatrix(presaved_path + "bact_filenames.txt")
else:
    print("Reading training data...")
    known_viruses, virus_filenames, virus_ids, virus_species = ReadDataBaseFilenames(virus_database_path, 0, presaved_path + "virus_filenames.txt", "virus_species.txt", virus_extension)
    known_bacterias, bact_filenames, bact_ids, bact_species = ReadDataBaseFilenames(bact_database_path, 0, presaved_path + "bact_filenames.txt", "bact_species.txt", bact_extension)
    virus_matrix = []
    bact_matrix = []
    n = 4; D = CreateDictionary(n)
    print("Generating viruses feature vectors...")
    for w in known_viruses:
        arr = FeatureVector(D, str(w), n)
        arr = np.divide(np.array(arr), len(w))
        virus_matrix.append(arr)
    m=0
    print("Generating bacterias feature vectors...")
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

matrix = np.vstack((virus_matrix,bact_matrix))
virus_bact_filenames = np.hstack((virus_filenames,bact_filenames))
virus_bact_species = np.hstack((virus_species,bact_species))

print("Reading lake files...")
lake, lake_filenames, lake_ids, lake_species = ReadDataBaseFilenames(lake_path, lake_quantity, presaved_path + "lake_filenames.txt", "lake_species.txt", lake_extension)
print("Read", len(lake), "file samples")
indices, dist = process()

saveOutput()
