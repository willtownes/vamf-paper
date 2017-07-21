"""
Find all FASTQ files in a specified directory (does not check subfolders). The fastq file names must have a certain structure. The suffix must be either ".fastq" or ".fastq.gz", case insensitive. The indication of paired end reads must be marked by an underscore followed by a substring containing "1" or "2". For example, a valid filename is "SRR123_r1.fastq" or "bc86_a01_2.fastq.gz".
Then, run kallisto on each set of paired end reads. Kallisto output is placed into a specified output folder with one subfolder for each of the fastq file pairs.
If a transcriptome file is not specified, the default is Homo sapiens GRCh38 (downloaded from Kallisto website), which the script searches for in the current directory and if not found it downloads and compiles the index.

Input: a folder containing a bunch of paired-end fastq files (first command line argument)
(optional) Output: a folder to store the output, where the original folder's structure will be replicated but will have kallisto output instead of fastq. By default the output folder is "kallisto_out" (second command line argument)

WARNING: Paired end reads only!
WARNING: If you have multiple files per sample, merge them first using fastq_merge.py before running this script.

Typical syntax
$ python ../util/kallisto_wrapper.py data/original/BC-84-P12/fastq_files_merged data/original/kallisto_out/BC-84-P12 -s "Homo sapiens"

original author: Will Townes (will.townes@gmail.com)
"""

#features to add in the future:
#(optional) location of transcriptome file used by kallisto
#(optional) is_paired_end is a true/false flag. Defaults to true.

import requests
import subprocess
import shlex
from os import path,listdir
from argparse import ArgumentParser
from misc import mkdir_p

CMD_BASE = "bsub -q medium -n 1 -R rusage[mem=2048] -o {bsub_out}/{fastq_id}.out -e {bsub_out}/{fastq_id}.err -J {fastq_id}"
#FOLDR=path.join("data/original")
#FOLDR=path.join("/data/aryee/ellisen/tnbc")
#CMD_BASE = "bsub -q medium -n 1 -R rusage[mem=2048] -o bsub_out/{sample}_{cell}.out -e bsub_out/{sample}_{cell}.err -J {sample}_{cell}"
# from http://bio.math.berkeley.edu/kallisto/transcriptomes/
SPECIES = {"Arabidopsis thaliana":"Arabidopsis_thaliana.TAIR10.26.cdna.all.fa.gz",
    "Caenorhabditis elegans":"Caenorhabditis_elegans.WBcel235.rel79.cdna.all.fa.gz",
    "Danio rerio":"Danio_rerio.Zv9.rel79.cdna.all.fa.gz",
    "Drosophila melanogaster":"Drosophila_melanogaster.BDGP6.rel79.cdna.all.fa.gz",
    "Homo sapiens":"Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz",
    "Latimeria chalumnae":"Latimeria_chalumnae.LatCha1.rel79.cdna.all.fa.gz",
    "Mus musculus":"Mus_musculus.GRCm38.rel79.cdna.all.fa.gz",
    "Rattus norvegicus":"Rattus_norvegicus.Rnor_5.0.rel79.cdna.all.fa.gz",
    "Saccharomyces cerevisiae":"Saccharomyces_cerevisiae.R64-1-1.rel81.cdna.all.fa.gz"}

class Kallisto_Wrapper_Error(Exception):
    """error raised if conditions required by this module not satisfied"""
    pass

def arghandler(args=None):
    '''parses a list of arguments (default is sys.argv[1:]), such as those from sys.argv and returns the parameters needed for quantifying fastq files with kallisto.'''
    helpmsg = '''This kallisto wrapper script takes as input a directory with FASTQ files and a species name, runs kallisto quantification on all files in parallel using the BSUB (LSF) system, then outputs the results into the specified output folder.'''
    parser = ArgumentParser(description=helpmsg)
    parser.add_argument("input_folder",type=str,help="Location of the folder containing fastq.gz files")
    parser.add_argument("output_folder",type=str,help="Location of the folder for storing Kallisto output")
    parser.add_argument("-s","--species",type=str,default="Homo sapiens",dest="species",help="Scientific name for the species of all the fastq files; used to indicate which transcriptome should be used. Scientific Names currently parsed include: {}. Default is 'Homo sapiens'".format(SPECIES.keys()))
    parser.add_argument("-t","--transcriptome-folder",type=str,default="../resources",dest="tf",help="Folder where transcriptomes are stored. Defaults to ../resources")
    parser.add_argument("-k","--kmer-size",type=int,default=None,dest="kmer_size",help="k-mer (odd) length, default is set to kallisto default. Must be less than fragment length")
    parser.add_argument("-n","--names",type=str,default=None,dest="sample_names",help="Path to file containing a list of sample names to be processed. If omitted, all samples in input folder are processed. Should have one sample name per row and no header")
    parser.add_argument('-v','--verbose',type=bool,default=False,dest="verbose",help="flag for printing detailed output")
    args = parser.parse_args(args) #if args is None, this will automatically parse sys.argv[1:]
    return args

#generate list of paths to fastq files
def is_fastq(filename):
    """returns true if the file is either a fastq file or a fastq.gz file"""
    split1 = path.splitext(filename)
    suffix = split1[-1].lower()
    if suffix==".fastq":
        return True
    elif suffix==".gz" and path.splitext(split1[0])[-1].lower()==".fastq":
        return True
    else:
        return False

def species2transcriptomeindex(species,kmer_size=None):
    """given a species name, return the name of the associated Kallisto index file"""
    spname = species.replace(" ","_")
    if kmer_size is None:
        return spname+".idx"
    else:
        return spname+"_kmer_"+str(kmer_size)+".idx"

def get_fastq_dict(foldr,paired_end_char="_"):
    """returns a dictionary of all fastq or fastq.gz files within the foldr specified. The keys are the identifiers of the fastq files. Each key has a list of length two associated with it, corresponding to the paired end reads. The "paired_end_char" indicates which character is used to separate the read1/read2 indicator information in the filename, just before the file extension. For example, if the fastq file name is SRR123_abc_r1.fastq.gz, then the paired end char is the underscore and the identifier of the file is SRR123_abc.
    The structure of the output is of the form {"sra_id":{"r1":"abc_1.fastq","r2":"abc_2.fastq"}}"""
    #note that we do not check to make sure there are exactly two reads per ID. New feature to be added later.
    res = {}
    for i in listdir(foldr):
        if not is_fastq(i): continue
        sra_id,tail = i.rsplit(paired_end_char,1)
        read = tail.split(".",1)[0]
        if "1" in read: flag="r1"
        elif "2" in read: flag="r2"
        else: raise Kallisto_Wrapper_Error("FASTQ Filename %s does not appear to have a paired end indicator 1 or 2"%i)
        try: res[sra_id][flag] = i
        except KeyError: res[sra_id] = {flag: i}
    return res

def check_transcriptome(folder,species,kmer_size=None):
    """checks in the specified folder for a valid kallisto transcriptome. If not found, downloads the latest version of the transcriptome corresponding to specified species to that folder"""
    print("Checking folder '%s' for existing kallisto transcriptome index file"%folder)
    transcriptome = SPECIES[species]
    tr_idx = path.join(folder,species2transcriptomeindex(species,kmer_size))
    #print(tr_idx)
    tr_gz = path.join(folder,transcriptome)
    if not path.exists(tr_idx):
        if not path.exists(tr_gz): #download transcriptome file if not exist
            mkdir_p(path.normpath(folder))
            print("%s not found, downloading..."%tr_gz)
            #wget http://bio.math.berkeley.edu/kallisto/transcriptomes/Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz
            r = requests.get("http://bio.math.berkeley.edu/kallisto/transcriptomes/"+transcriptome,stream=True)
            with open(tr_gz,"wb") as fd:
                for chunk in r.iter_content():
                    fd.write(chunk)
        # compile kallisto index if not exist
        print("%s not found, regenerating..."%tr_idx)
        #kallisto index -i transcripts.idx Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz
        if kmer_size is None: kallisto_args = ["kallisto","index","-i",tr_idx,tr_gz]
        else: kallisto_args = ["kallisto","index","-i",tr_idx,"-k",str(kmer_size),tr_gz]
        subprocess.call(kallisto_args)
    return tr_idx

def quant(input_folder,fastq_dict,output_folder,transcriptome,bsub_out_toplevel="bsub_out"):
    """submit jobs to bsub to run kallisto quantification on each pair of fastq files in the fastq dictionary provided, where the files are located in the input_folder."""
    print("Starting new quantification run for batch of %d read pairs from %s"%(len(fastq_dict),input_folder))
    mkdir_p(bsub_out_toplevel)
    #hack to recover batch name, need to update to obtain from function arguments.
    batch_name = path.split(output_folder)[-1]
    print("Inferred batch name: %s"%batch_name)
    bsub_out = path.join(bsub_out_toplevel,batch_name)
    mkdir_p(bsub_out)
    print("bsub logs stored in %s folder"%bsub_out)
    mkdir_p(output_folder)
    print("kallisto output in %s"%output_folder)
    for i in fastq_dict:
        print("===processing fastq files from FASTQ ID: %s==="%i)
        outdir = path.join(output_folder,i) #separate folder for each fastq, within the output folder
        mkdir_p(outdir)
        cmd = CMD_BASE.format(fastq_id=i,bsub_out=bsub_out)
        cmd = shlex.split(cmd) #convert to list of arguments
        try:
            f1 = path.join(input_folder,fastq_dict[i]["r1"])
            f2 = path.join(input_folder,fastq_dict[i]["r2"])
        except KeyError:
            print("Error- missing paired end read file!")
        else:
            cmd.append("kallisto quant -i {index} -o {out} {fastq_1} {fastq_2}".format(index=transcriptome,out=outdir,fastq_1=f1,fastq_2=f2))
            subprocess.call(cmd)

if __name__=="__main__":
    #command line parsing
    args = arghandler()
    tr_idx = check_transcriptome(args.tf,args.species,args.kmer_size)
    fastq_dict = get_fastq_dict(args.input_folder)
    if args.sample_names is not None:
      sample_names = [i.strip() for i in open(args.sample_names).readlines()]
      try:
        fastq_dict = dict((k,fastq_dict[k]) for k in sample_names)
      except KeyError:
        raise Kallisto_Wrapper_Error("Sample names must be a subset of files in input folder")
    quant(args.input_folder,fastq_dict,args.output_folder,tr_idx)
