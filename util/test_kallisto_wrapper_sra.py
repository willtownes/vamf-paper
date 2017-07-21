"""
Tests for kallisto_wrapper_sra module
Recommended command line syntax:
py.test --cov=.
"""

import kallisto_wrapper_sra as kws

fqdict = {"s1":{"fastq_list":["srr123","srr456"],"is_paired_end":True,"species":"Homo sapiens","avgLengths":[50,80]},
            "s2":{"fastq_list":["srr789"],"is_paired_end":False,"species":"Homo sapiens","avgLengths":[30]},
            "s3":{"fastq_list":["srr333"],"is_paired_end":True,"species":"Mus musculus","avgLengths":[150]}}

def test_1species_mixed_library():
    fqd = {i: fqdict[i] for i in fqdict if fqdict[i]["species"] == "Homo sapiens"}
    assert kws.compute_ideal_kmer(fqd) == {"Homo sapiens":11}

def test_2species_same_library():
    fqd = {i: fqdict[i] for i in fqdict if i in ["s1","s3"]}
    assert kws.compute_ideal_kmer(fqd) == {"Homo sapiens":11,"Mus musculus":None}

def test_single_end():
    fqd = {i: fqdict[i] for i in fqdict if i == "s2"}
    assert kws.compute_ideal_kmer(fqd) == {"Homo sapiens":15}
