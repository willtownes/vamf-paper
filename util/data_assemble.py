"""
Read sample and cell-related QC information from bsub_out and save in a summary file.
Joins the pseudo-alignments for the different cells into a single combined abundance file (csv)
INPUT PARAMETERS:
* bsub_folder=bsub_out
* kallisto_folder=data/original/kallisto_out
* sample=None
* output_folder=data
OUTPUT:
* two output files stored in a folder named [output_folder]/[sample]
* Concatenated alignments are placed in a file called "abundance_combined.csv"
* combined QC results are placed in a file called "kallisto_qc.csv"

Command line usage:
python ../util/data_assemble.py -s BC1PT39
"""
import os,csv,re
from argparse import ArgumentParser
from misc import mkdir_p

class Data_Assemble_Error(Exception):
    """error raised if conditions required by this module not satisfied"""
    pass

def get_readcounts(filepath):
    """takes a kallisto output file string (in bsub_out folder, a file like bc1_A01.err), and extracts the total reads and aligned reads, returned as a tuple"""
    for line in open(filepath):
        if "[quant] processed " and "reads pseudoaligned" in line:
            l = line.split("[quant] processed ")[1].split(" reads pseudoaligned")[0]
            total,aligned = l.split(" reads, ")
            total = int("".join(total.split(",")))
            aligned = int("".join(aligned.split(",")))
            return total,aligned #terminate once the line is found

def lstrip_batchname(batch_name,cell_name):
    """strip out the extra batch name sometimes present in beginning of cell names for more concise notation"""
    if batch_name=="":
      return cell_name
    else:
      pattern = re.compile(batch_name+'[-_]*') #batch name, optionally followed by - or _
      return pattern.split(cell_name)[-1]

def get_kallisto_qc(batch_name,bsub_dir,verbose=True):
    """reads the Kallisto logs from the bsub output folder, parses QC information, returns as list of tuples suitable for writing to CSV. Must specify the batch name (ie, the cluster of cells which were processed together, often from a single tumor)"""
    #bsub_dir = os.path.join(bsub_dir,batch_name)
    flist = [i for i in os.listdir(bsub_dir) if batch_name in i]
    outs = [i for i in flist if ".out" in i]
    errs = [i for i in flist if ".err" in i]
    if len(outs)!=len(errs):
      raise Data_Assemble_Error("BSUB .out files don't match with BSUB .err files!")
    res = [("sample","total_reads","aligned_reads")]
    for i in xrange(len(outs)):
        status = open(os.path.join(bsub_dir,outs[i])).read()
        sample = outs[i].split(".")[0]
        #remove redundant batch label from sample id, if exists
        #sample = lstrip_batchname(batch_name,sample)
        if "Successfully completed." not in status:
            if verbose: print("Sample %s failed to process"%(outs[i].split(".")[0]))
            continue
        else:
            total,aligned = get_readcounts(os.path.join(bsub_dir,errs[i]))
            res.append((sample,total,aligned))
    return res

def write_kallisto_qc_to_csv(qc_data,csv_file):
    """write out list of tuples to the specified file. Will overwrite file if already existing"""
    with open(csv_file,"w") as ofile:
        writer = csv.writer(ofile)
        for row in qc_data:
            writer.writerow(row)

def stack_counts(batch_name,indir="data/original/kallisto_out",outdir="data/",verbose=True,eps=.0001):
    """read each kallisto output file from a batch and combine them into a single file,
    omitting all the zero counts"""
    kallisto_header = ['target_id', 'length', 'eff_length', 'est_counts', 'tpm']
    header = ['cell']+kallisto_header
    with open(os.path.join(outdir,batch_name,"abundance_combined.csv"),"w") as ofile:
        writer = csv.DictWriter(ofile,header)
        writer.writeheader()
        for i in os.listdir(indir): #i is the name of a cell
            if batch_name not in i: continue
            #cell_label = {'cell':lstrip_batchname(batch_name,i)}
            cell_label = {'cell':os.path.splitext(i)[0]}
            if verbose: print("Now stacking abundances from %s to combined file"%cell_label['cell'])
            with open(os.path.join(indir,i,"abundance.tsv"),"r") as ifile:
                reader = csv.DictReader(ifile,dialect="excel-tab")
                for line in reader:
                    if float(line["est_counts"])<eps: continue
                    else:
                        line.update(cell_label)
                        writer.writerow(line)

def arghandler(args=None):
    '''parses a list of arguments (default is sys.argv[1:]), such as those from sys.argv and returns the parameters needed for aligning fastq files with kallisto.'''
    helpmsg = '''Script to combine results from multiple cells in the same sample. Creates two files: abundance_combined.csv and kallisto_qc.csv'''
    parser = ArgumentParser(description=helpmsg)
    parser.add_argument("-b","--bsub-folder",type=str,default="bsub_out",dest="bsub",help="Location of the folder containing output from LSF bsub system, to extract quality control information. Default: 'bsub_out'")
    parser.add_argument("-k","--kallisto-folder",type=str,default=os.path.join('data','original','kallisto_out'),dest="kallisto",help="Location of the folder containing Kallisto output. Default: 'data/original/kallisto_out'")
    parser.add_argument("-s","--sample-name",type=str,default="",dest="sample",help="Optional sample name to combine only results with names containing this as a substring. Default: combine all results")
    parser.add_argument("-o","--output-folder",type=str,default="data",dest="output_folder",help="Folder where combined results are stored")
    parser.add_argument('-v','--verbose',type=bool,default=True,dest="verbose",help="flag for printing detailed output")
    args = parser.parse_args(args) #if args is None, this will automatically parse sys.argv[1:]
    return args
    
if __name__=="__main__":
    args = arghandler()
    out = os.path.join(args.output_folder,args.sample)
    mkdir_p(out)
    qc_dat = get_kallisto_qc(args.sample,args.bsub,args.verbose)
    write_kallisto_qc_to_csv(qc_dat,os.path.join(out,"kallisto_qc.csv"))
    stack_counts(args.sample,args.kallisto,args.output_folder,args.verbose)




