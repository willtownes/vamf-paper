"""
Run fSCLVM on data and write out the estimated sparse factors for the cells.
Input data should have genes in the rows, cells/samples in the columns (note this is the transpose of the format specified in the fsclvm documentation, but is more standard for gene expression data).
It is assumed that the data has already been filtered to remove genes that are expressed in <20 percent of samples.
Also, it is assumed that the data have already been normalized (and usually it's good to convert to log scale as well)
Optionally, the data can have a header row, which will be removed before running fsclvm.
The script writes out the low-dimensional projection of the cells' expression profiles, as well as the model parameters.
fSCLVM can be installed using pip install fsclvm

example command
time python fsclvm_wrapper.py data/noise_only.tsv --log-transform False --nHiddenSparse 2
"""

import fscLVM
import numpy as np
from scipy import sparse as sps
from os import path
from argparse import ArgumentParser

#class FSCLVM_Error(Exception):
#    pass

def read_data(filename, logtrans=True, delim=" ",skiprows=1):
    """Read in the coo data file at filename and return a sparse array suitable for processing in fsclvm.
    CSV has first column row ID, second column col ID, third column data value.
    No rownames allowed. Checks if zero or one based indexing is used."""
    x = np.loadtxt(filename,delimiter=delim,skiprows=skiprows)
    i = x[:,0]
    j = x[:,1]
    #adjust for non-zero indexing
    i = i-min(i)
    j = j-min(j)
    if(logtrans): x[:,2] = np.log2(1+x[:,2])
    return sps.coo_matrix( (x[:,2], (i,j)) ).tocsc()

def fsclvm_wrapper(X,**kwargs):
    """run fsclvm on the numpy array in data. Assumes that the data have been normalized.
    X should have samples (cells) in the columns and features (genes) in the rows.
    **kwargs passed to fSCLVM initFA function (eg, nHidden, nHiddenSparse, noise, etc)"""
    I = np.ones((X.shape[0],1)) #length= number of genes
    terms = np.array(["dense0"],dtype='|S21') #one fake pathway
    FA = fscLVM.initFA(X.T.toarray(),terms,I,**kwargs)
    FA.train()
    #return FA.getX()
    return FA

def extract_factors(FA,sparse=True):
    """extract either the sparse or the dense factors from factor analysis object produced by fsclvm_wrapper function.
    If sparse==True, returns the sparse factors, otherwise, returns the dense factors."""
    terms = FA.getTerms()
    if sparse: terms = sorted([x for x in terms if 'sparse' in x.lower()])
    else: terms = sorted([x for x in terms if 'sparse' not in x.lower()])
    return FA.getX(terms=terms)

def arghandler(args=None):
    '''parses a list of arguments (default is sys.argv[1:])'''
    helpmsg = '''This fSCLVM wrapper script takes as input a tab-delimited file with sparse matrix data (row,column,data). The first row is assumed to be a header and ignored. fSCLVM is then run on the matrix and the factors are written to disk.'''
    parser = ArgumentParser(description=helpmsg)
    parser.add_argument("ifile",type=str,help="Location of the input file containing the sparse matrix")
    parser.add_argument("-o","--ofile",type=str,default=None,dest="ofile",help="Location of the file for storing ZIFA output (the factors). Default is same file name as input file with '_fsclvm.tsv' appended")
    parser.add_argument("-k","--nHidden",type=int,default=2,dest="nHidden",help="Desired dimensionality of dense factors, defaults to 2")
    parser.add_argument("-s","--nHiddenSparse",type=int,default=0,dest="nHiddenSparse",help="Desired dimensionality of sparse factors, defaults to 2")
    parser.add_argument("-n","--noise",type=str,default="gauss",dest="noise",help="fSCLVM noise model, can be 'gauss', 'hurdle', or 'poisson'")
    parser.add_argument("-l","--log-transform",type=bool,default=True,dest="logtrans",help="If True, the data values are transformed by log2(1+x) before applying the fSCLVM method")
    parser.add_argument("-d","--dense-save",type=bool,default=False,dest="dsave",help="Should dense factors be saved to disk? Default: saves the sparse factors only")
    args = parser.parse_args(args) #if args is None, this will automatically parse sys.argv[1:]
    if args.ofile is None:
        args.ofile = path.splitext(args.ifile)[0]+"_fsclvm.tsv"
    return args

if __name__=="__main__":
    args = arghandler()
    dat = read_data(args.ifile, logtrans=args.logtrans) #sparse csc format
    FA = fsclvm_wrapper(dat, nHidden=args.nHidden, nHiddenSparse=args.nHiddenSparse, noise=args.noise)
    factors = extract_factors(FA,sparse=(not args.dsave))
    np.savetxt(args.ofile,factors)

    #dat = read_data("data/latent_clusters.tsv",logtrans=False)
    #FA = fsclvm_wrapper(dat, nHidden=2, nHiddenSparse=2, noise="gauss")
    #factors_sparse = extract_factors(FA,sparse=True)
    #factors_dense = extract_factors(FA,sparse=False)
    #from matplotlib import pyplot as plt
    #plt.scatter(factors_sparse[:,0],factors_sparse[:,1])
    #plt.scatter(factors_dense[:,0],factors_dense[:,2])


