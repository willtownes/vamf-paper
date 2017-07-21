"""
You have a bunch of GEO identifier named subfolders of kallisto_output after quantification of a bunch of FASTQ files. 
You want to merge the kallisto outputs based on which cells are in the same grouping (eg tumor).
"""

import csv
from os import path,symlink
from argparse import ArgumentParser
from misc import mkdir_p

pj = path.join
def arghandler(args=None):
    """parse command line arguments"""
    helpmsg="create symbolic links to include group names as well as GEO IDs in output folder of an RNA-Seq quantifier such as Kallisto"
    parser = ArgumentParser(description=helpmsg)
    parser.add_argument(type=str,dest="mapfile",help="file path to a CSV file (assumed to have a header in first row). CSV column 1 contains the GEO IDs (unique sample IDs). CSV column 2 contains the corresponding group name (eg tumor)")
    parser.add_argument(type=str,dest="quant_folder",help="directory path to the parent folder of the GEO identified quantification output folders (eg kallisto_out)")
    parser.add_argument("-b","--bsub_folder",type=str,default="bsub_out",dest="bsub_folder",help="directory path to the bsub output folder, containing logfiles of the quantifier")
    #parser.add_argument("-a","--alias_folder",type=str,default=None,dest="alias_folder",help="directory path to where the symbolic links should be stored. Default is to add '_alias' to quant_folder")
    args = parser.parse_args(args)
    return args

if __name__=="__main__":
    args = arghandler()
    #if args.alias_folder is None:
    alias_folder = args.quant_folder + "_alias"
    bsub_alias = args.bsub_folder + "_alias"
    mkdir_p(alias_folder)
    mkdir_p(bsub_alias)
    with open(args.mapfile) as f:
        next(f) #skip header
        reader = csv.reader(f)
        for line in reader:
            s = line[0]
            ln_name = line[1]+"_"+s
            symlink(path.abspath(pj(args.quant_folder,s)), pj(alias_folder,ln_name))
            symlink(path.abspath(pj(args.bsub_folder,s+".err")), pj(bsub_alias,ln_name+".err"))
            symlink(path.abspath(pj(args.bsub_folder,s+".out")), pj(bsub_alias,ln_name+".out"))
