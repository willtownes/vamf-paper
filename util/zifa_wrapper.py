"""
Run ZIFA on data and write out the estimated factors for the cells.
Input data should have genes in the rows, cells/samples in the columns (note this is the transpose of the format specified in the ZIFA documentation, but is more standard for gene expression data).
It is assumed that the data has already been filtered to remove genes that are expressed in <20 percent of samples.
Also, it is assumed that the data have already been normalized (and usually it's good to convert to log scale as well)
Optionally, the data can have a header row and/or the first column containing gene IDs. These will be removed before running ZIFA.
The script writes out a pickled representation of the low-dimensional projection of the cells' expression profiles, as well as the model parameters.
ZIFA can be installed using pip install git+https://github.com/epierson9/ZIFA.git#egg=ZIFA
"""

from ZIFA import ZIFA,block_ZIFA
import numpy as np
from scipy import sparse as sps
from sys import argv
from os import path
from argparse import ArgumentParser

def read_data(filename, logtrans=True, delim=" ",skiprows=1):
    """Read in the coo data file at filename and return a sparse array suitable for processing in ZIFA. CSV has first column row ID, second column col ID, third column data value. No rownames allowed. Checks if zero or one based indexing is used."""
    x = np.loadtxt(filename,delimiter=delim,skiprows=skiprows)
    x[:,2][x[:,2]<0] = 0 #handle case of negative values
    i = x[:,0]
    j = x[:,1]
    #adjust for non-zero indexing
    i = i-min(i)
    j = j-min(j)
    if(logtrans): x[:,2] = np.log2(1+x[:,2])
    return sps.coo_matrix( (x[:,2], (i,j)) ).tocsc()

def zifa_wrapper(X, k=2, block=True, return_model_params=False, p0_thresh=None):
    """run ZIFA on the numpy array in data. Assumes that the data have been normalized. X should have samples (cells) in the columns and features (genes) in the rows."""
    if block:
        if p0_thresh is None: Z, model_params = block_ZIFA.fitModel(X.T.toarray(), k)
        else: Z, model_params = block_ZIFA.fitModel(X.T.toarray(), k)
    else:
        Z, model_params = ZIFA.fitModel(X.T.toarray(), k)
    if return_model_params: return Z,model_params
    else: return Z

#def write_data(Z, ofile, model_params = None):
#    """write out the factors and optionally model parameters for later use.
#    Output file names will have ofile_base followed by _z.tsv for Z and _params_pickle.txt for model_params"""
    #np.savetxt(ofile_base+"_zifa.tsv",Z)
    #if model_params is not None:
    #    with open(ofile_base+"_zifa_params_pickle.txt","w") as ofile:
    #        pickle.dump(model_params,ofile)

def arghandler(args=None):
    '''parses a list of arguments (default is sys.argv[1:]), such as those from sys.argv and returns the parameters needed for running ZIFA on a '''
    helpmsg = '''This ZIFA wrapper script takes as input a tab-delimited file with sparse matrix data (row,column,data). The first row is assumed to be a header and ignored. ZIFA is then run on the matrix and the factors are written to disk.'''
    parser = ArgumentParser(description=helpmsg)
    parser.add_argument("ifile",type=str,help="Location of the input file containing the sparse matrix")
    parser.add_argument("-o","--ofile",type=str,default=None,dest="ofile",help="Location of the file for storing ZIFA output (the factors). Default is same file name as input file with '_z.tsv' appended")
    parser.add_argument("-k","--dim",type=int,default=2,dest="k",help="Desired dimensionality of factors, defaults to 2 for visualization")
    parser.add_argument("-b","--block",type=bool,default=False,dest="block",help="Flag for using the block ZIFA approximate algorithm rather than the full ZIFA. Block ZIFA is faster but less accurate.")
    parser.add_argument("-p","--p0-thresh",type=float,default=None,dest="p0_thresh",help="If block ZIFA is used, what threshold to use for filtering genes with many zeroes. ZIFA automatically sets this to 0.95 if not provided. Recommend not to change unless ZIFA throws a warning or has convergence problems")
    parser.add_argument("-l","--log-transform",type=bool,default=True,dest="logtrans",help="If True, the data values are transformed by log2(1+x) before applying the ZIFA method")
    args = parser.parse_args(args) #if args is None, this will automatically parse sys.argv[1:]
    if args.ofile is None:
        args.ofile = path.splitext(args.ifile)[0]+"_zifa.tsv"
    return args

if __name__=="__main__":
    args = arghandler()
    dat = read_data(args.ifile, logtrans=args.logtrans) #sparse csc format
    z = zifa_wrapper(dat, k=args.k, block=args.block, p0_thresh=args.p0_thresh)
    np.savetxt(args.ofile,z)
