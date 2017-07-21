"""
Tests for fastq_merge module
Recommended command line syntax:
py.test --cov=.
"""

import fastq_merge as fm

def test_filename2meta_1(): #normal case
    test_str = "BC-84-P12-D10_S291_L002_R1_002.fastq.gz"
    correct = {"sample":"BC-84-P12-D10_S291","lane":1,"read":"r1","file":1,"gz":True}
    assert fm.filename2meta(test_str)==correct

def test_filename2meta_2(): #no file id
    test_str = "BC-84-P12-D10_S291_L001_R1.fastq.gz"
    correct = {"sample":"BC-84-P12-D10_S291","lane":0,"read":"r1","file":0,"gz":True}
    assert fm.filename2meta(test_str)==correct

def test_filename2meta_3(): #no lane id
    test_str = "BC-84-P12-D10_S291_R1_001.fastq.gz"
    correct = {"sample":"BC-84-P12-D10_S291","lane":0,"read":"r1","file":0,"gz":True}
    assert fm.filename2meta(test_str)==correct

def test_filename2meta_4(): #read 2, lower case
    test_str = "BC-84-P12-D10_S291_L001_r2_001.fastq.gz"
    correct = {"sample":"BC-84-P12-D10_S291","lane":0,"read":"r2","file":0,"gz":True}
    assert fm.filename2meta(test_str)==correct

def test_filename2meta_5(): #non gz file
    test_str = "BC-84-P12-D10_S291_L001_r2_001.fastq"
    correct = {"sample":"BC-84-P12-D10_S291","lane":0,"read":"r2","file":0,"gz":False}
    assert fm.filename2meta(test_str)==correct

def test_filename2meta_6(): #alternative name scheme
    test_str = "F09_R1.fastq"
    correct = {"sample":"F09","lane":0,"read":"r1","file":0,"gz":False}
    assert fm.filename2meta(test_str)==correct

#need to add tests for the cases that throw exceptions

file_list=[
    "BC-84-P12-D11_S292_L001_R2_001.fastq.gz",
    "BC-84-P12-D10_S291_L001_R1_002.fastq.gz",
    "BC-84-P12-D10_S291_L001_R1_001.fastq.gz",
    "BC-84-P12-D11_S292_L001_R1_001.fastq.gz",
    "BC-84-P12-D10_S291_L001_R2_001.fastq.gz",
    "BC-84-P12-D10_S291_L002_R1_001.fastq.gz",
    "BC-84-P12-D10_S291_L002_R2_001.fastq.gz"] #note files are out of order


# expected_file_dict={
# "BC-84-P12-D10_S291":{
#     "r1":{
#         (0,0):"BC-84-P12-D10_S291_L001_R1_001.fastq.gz",
#         (0,1):"BC-84-P12-D10_S291_L001_R1_002.fastq.gz",
#         (1,0):"BC-84-P12-D10_S291_L002_R1_001.fastq.gz"},
#     "r2":{(0,0):"BC-84-P12-D10_S291_L001_R2_001.fastq.gz",
#         (1,0):"BC-84-P12-D10_S291_L002_R2_001.fastq.gz"}
#     },
# "BC-84-P12-D11_S292":{
#     "r1":{(0,0):"BC-84-P12-D11_S292_L001_R1_001.fastq.gz"},
#     "r2":{(0,0):"BC-84-P12-D11_S292_L001_R2_001.fastq.gz"}}
# }

expected_file_dict={
"BC-84-P12-D10_S291":{
    "r1":("BC-84-P12-D10_S291_L001_R1_001.fastq.gz",
        "BC-84-P12-D10_S291_L001_R1_002.fastq.gz",
        "BC-84-P12-D10_S291_L002_R1_001.fastq.gz"),
    "r2":("BC-84-P12-D10_S291_L001_R2_001.fastq.gz",
        "BC-84-P12-D10_S291_L002_R2_001.fastq.gz")
    },
"BC-84-P12-D11_S292":{
    "r1":("BC-84-P12-D11_S292_L001_R1_001.fastq.gz",),
    "r2":("BC-84-P12-D11_S292_L001_R2_001.fastq.gz",)}
}

def test_dict2tuple():
    assert fm.dict2tuple({(0,11):"abc",(0,1):"def"})==("def","abc")

def test_build_fastq_dict_1():
    assert fm.build_fastq_dict(file_list)==expected_file_dict

#def test_sort_verify_fastq_dict():
#    assert fm.sort_verify_fastq_dict(expected_file_dict)==expected_file_dict2
