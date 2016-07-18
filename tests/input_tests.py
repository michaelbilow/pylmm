from pylmm.input import plink
from pylmm.lmm import LMM
from os.path import join, split, splitext
import numpy as np
import sys

def test1(plink_object):
    ## Ensure the number of individuals == the number of lines in the FAM file
    counter = 0
    if plink_object.type == 'b':
        fam_file = plink_object.fbase + '.fam'
    elif plink_object.type == 't':
        fam_file = plink_object.fbase + '.tfam'
    with open(fam_file, 'r') as opfile:
        for line in opfile:
            counter += 1
    assert counter == len(plink_object.indivs), \
        "getIndivs failure: counted {} lines, got  {} indivs".format(counter, len(plink_object.indivs))
    return

def test2(plink_object):
    ## Ensure the number of individuals == the number of phenotypes read in
    assert plink_object.phenos is not None, "Please specify a phenotype file."
    assert len(plink_object.phenos) == len(plink_object.indivs), \
        "Number of phenotypes {} does not equal number of individuals read {}".format(plink_object.phenos, plink_object.indivs)
    return

def test3(plink_object):
    ## Check that a kinship file has been read
    assert plink_object.kFile is not None, "Please specify a kinship matrix"
    return

def test4(plink_object):
    ## Check that each individuals phenotype corresponds correctly with the phenofile
    indiv_index = 0
    with open(plink_object.phenoFile, 'r') as opfile:
        for line in opfile:
            assert round(float(line.strip().split()[2]), 10) == round(plink_object.phenos[indiv_index][0], 10), \
                "Phenotype read incorrectly for individual {}".format(indiv_index)
            indiv_index += 1
    return

def test5(plink_object):
    ## Check SNPs in plink object are the same as those in the base file
    if plink_object.type != 't':
        return
    with open(plink_object.fbase + '.tped', 'r') as f:
        for snp, id in plink_object:
            line = f.readline().strip().split()[4:]
            line_pairs = zip(line[::2], line[1::2])
            for snp_val, line_pair in zip(snp, line_pairs):
                if all(x == '2' for x in line_pair):
                    assert snp_val == 1.0, "SNPs read incorrectly {}".format(line)
                elif all(x == '1' for x in line_pair):
                    assert snp_val == 0.0, "SNPs read incorrectly {}".format(line)
                elif all(x == '0' for x in line_pair):
                    assert np.isnan(snp_val), "SNPs read incorrectly {}".format(line)
                else:
                    assert snp_val == 0.5, "SNPs read incorrectly {}".format(line)
    return

def test6(plink_object):
    ## Check that the kinship matrix is read properly
    line_index = 0
    with open(plink_object.kFile, 'r') as opfile:
        for line in opfile:
            itr = line.strip().split()
            for i in range(len(plink_object.K[line_index,])):
                assert float(itr[i]) == plink_object.K[line_index, i], \
                    "Kinship matrix was not read properly"
            line_index += 1
    return

def run_all_tests(plink_object):
    test1(plink_object)
    test2(plink_object)
    test3(plink_object)
    test4(plink_object)
    test5(plink_object)
    test6(plink_object)
    return

if __name__ == "__main__":
    input_folder = '../data/'
    mouse_data = False
    if mouse_data:
        base_name = 'snps.132k.clean.noX{}'
        my_plink_object = plink(fbase=join(input_folder, base_name.format('')),
                                kFile=join(input_folder, base_name.format('.pylmm.kin')),
                                phenoFile=join(input_folder, base_name.format('.fake.phenos')), type='b',
                                normGenotype=True, readKFile=False, fastLMM_kinship=False)
    else:
    ### for finland human data, need pheno file though
        base_name = 'finland_short{}'
        my_plink_object = plink(fbase=join(input_folder, base_name.format('')),
                                kFile=join(input_folder, base_name.format('.pylmm.kin')),
                                phenoFile=None, type='t', normGenotype=False, readKFile=False,fastLMM_kinship=False)

    run_all_tests(my_plink_object)
    lmm_object = LMM(my_plink_object.phenos, my_plink_object.K, verbose=True)
