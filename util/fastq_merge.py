"""
Scan a folder containing fastq or fastq.gz files,
creates a nested dictionary of structure
{sample_id:{read1:[list of fastq for read 1],read2:[list of fastq for read 2]}}
This facilitates merging the fastq files into one fastq per end of each sample (two fastq files per sample
"""
import re,os,shutil
#from copy import deepcopy
from sys import argv

class Fastq_Merge_Error(Exception):
    """Base class for exceptions in this module"""
    pass

def filename2meta(f):
    """takes a fastq filename and returns a dictionary describing its metadata"""
    fsplit=os.path.splitext(f)
    if fsplit[-1].lower()==".gz":
        gz=True
        fsplit = os.path.splitext(fsplit[0])
    else:
        gz = False
    if fsplit[-1].lower()==".fastq":
        ftxt = fsplit[0].split("_") #[sample-id,lane,read-id,file-id] usually
    else:
        raise Fastq_Merge_Error("Invalid filename, must be either .fastq or .fastq.gz")
    file_id_str = re.match(re.compile(r"\d+"),ftxt[-1])
    if file_id_str: #case of fastq split across multiple files
        file_id = int(file_id_str.group())-1 #shift to zero index for python
        ftxt.pop() #remove file id from ftxt
    else: #case where only one file per fastq
        file_id = 0
    read_id = ftxt.pop().lower()
    if read_id not in ("r1","r2"):
        raise Fastq_Merge_Error("Invalid or missing paired end read ID, should be r1,r2,R1,or R2")
    lane_id_str = re.match(re.compile(r"l\d+"),ftxt[-1].lower())
    if lane_id_str:
        lane_id = int(ftxt[-1][1:])-1
        ftxt.pop()
    else:
        lane_id = 0
    sample_id = "_".join(ftxt)
    return {"sample":sample_id,"lane":lane_id,"read":read_id,"file":file_id,"gz":gz}

#sort key is by (lane,file) within each sample+read combination

def dict2tuple(d):
    """converts a dictionary like this {(0,0):"abc",(0,1):"def"} into tuple like this ("abc","def") following the sort order of the keys"""
    return tuple(d[i] for i in sorted(d))

def build_fastq_dict(file_list):
    """given a list of files, build a nested dictionary to show which files belong to which sample."""
    #if len(file_list)%2 != 0:
    #    raise Fastq_Merge_Error("Odd Number of Files- Probably one or more files are missing!")
    fd = {}
    if ".fastq.gz" in file_list[0].lower(): gz=True
    else: gz=False
    for i in file_list:
        m = filename2meta(i)
        if m["gz"] != gz:
            raise Fastq_Merge_Error("Inconsistent gzip pattern. Either convert all files to .fastq.gz or convert all to .fastq")
        sample = m["sample"]
        rd = m["read"]
        sortkey = (m["lane"],m["file"]) #tuple
        if sample not in fd:
            fd[sample] = {rd:{sortkey:i}}
        elif rd not in fd[sample]:
            fd[sample][rd] = {sortkey:i}
        elif sortkey not in fd[sample][rd]:
            fd[sample][rd][sortkey] = i
        else:
            raise Fastq_Merge_Error("Duplicate Lane/File ID combination for Fastq %s"%i)
    #converts the tuple keyed inner dictionary into a tuple with proper ordering.
    #also checks that there is at least one file for R1 and R2 for all samples
    for s in fd:
        #assert len(fd[s]["r1"])==len(fd[s]["r2"])
        assert len(fd[s]["r1"])>0
        fd[s]["r1"] = dict2tuple(fd[s]["r1"])
        fd[s]["r2"] = dict2tuple(fd[s]["r2"])
    return fd

# def sort_verify_fastq_dict(fd):
#     """Takes a nested dictionary output from build_fastq_dict() and
#     """
#     #does not check that same number of files in R1 as R2 for all samples
#     fd = deepcopy(fd)
#     for s in fd:
#         #assert len(fd[s]["r1"])==len(fd[s]["r2"])
#         assert len(fd[s]["r1"])>0
#         fd[s]["r1"] = dict2tuple(fd[s]["r1"])
#         fd[s]["r2"] = dict2tuple(fd[s]["r2"])
#     return fd

# script to remove unwanted substring from each file in a directory, merge multiple files into one fastq per sample+read

def merge_fastq(infolder,outfolder):
    """read in all fastq files in 'infolder', parse filenames into dictionary, then concatenate files from the same sample+read into single files in the 'outfolder' directory. It will not overwrite an existing outfolder"""
    files = [i for i in os.listdir(infolder) if ".fastq" in i.lower()]
    if ".fastq.gz" in files[0].lower(): gz=True
    else: gz=False
    fd = build_fastq_dict(files)
    os.mkdir(outfolder)
    for s in fd: #iterate over samples
        for r in ("r1","r2"): #two reads per sample
            ofile_name = s+"_"+r+".fastq.gz"
            with open(os.path.join(outfolder,ofile_name),"wb") as ofile:
                for ifile_name in fd[s][r]:
                    with open(os.path.join(infolder,ifile_name),"rb") as ifile:
                        shutil.copyfileobj(ifile,ofile) #appends to ofile

if __name__=="__main__":
    infolder = argv[1]
    outfolder = argv[2]
    merge_fastq(infolder,outfolder)