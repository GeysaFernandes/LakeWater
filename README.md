# A tool for optimizing the process of sequence analysis

## Motivation
Lake Michigan is very important for Chicago area and nearby cities, and more than that it's one of the largest bodies of freshwater on the planet. Dr. Catherine Putonti's lab collected samples of the lake Michigan's water and extracted DNA fragments of the microorganisms that populated the samples. Then Dr. Putonti began to try to identify those microorganisms by searching for the extracted DNA sequences in a database of known microorganisms. But analyzing and classifying those DNA sequences can be a very hard task due to fragmented DNA, mutations, extensive database and computationally expensive algorithms used to match DNA sequences. As an example Dr. Putonti's lab took three months to analyse nine lake samples, and she has more than a million. 

## Objective
One of the most popular tools to match sequences is called string alignment, a method that tries to find similar regions between two sequences. Since it is a very computationally expensive method we want to provide a reduced number of candidates, that are more likely to score a sequence match, to use string alignment. So rather that using alignment in every single sequence in the database, the researchers can use string alignment just on the most likely to match sequences, reducing the time to process the data.
We measured our accuracy based on [Dr. Catherine Putonti sample collection from Lake Michigan] (http://www.publish.csiro.au/paper/MF15172.htm).

## Methods
###Main methods

1. k-d Tree

    • is the standard method to find the nearest neighbors.

###Supportive methods

1. Feature Vector 

    • The feature vector is a score vector for all possible combinations of the nucleobases arranged in groups of k bases. Advised by Dr. Putonti we arrange the nucleotides by 3, 4 or 6.

2. PCA - Principal Component Analysis

    • We used PCA as a dimensionality reduction method

## Files
### lake_water_accuracy_test.ipynb
This is the script that tests out method. Since we had only real data without groundtruth, we created a simulated case test.
For our test, we extracted subsequences of known viruses and bacterias and considered them as lakes. A right answer of our code is when the original virus or bacteria sequence is among the candidates.

#### How to run
For this one, you have to have:

1. [Python 3.5](https://www.python.org/downloads/)
2. [Jupyter notebook](http://jupyter.readthedocs.io/en/latest/install.html)
3. [BioPython](http://biopython.org/wiki/Download)
4. [SciKit Learn](http://scikit-learn.org/stable/)

### alignment_candidates.py
For this one, you will the same things from above but Jupyter notebook.
The parameters for the code are specified in the code.
To run, you just have to type the following snipped on Terminal/prompt:
```
$ python alignment_candidates.py
```
