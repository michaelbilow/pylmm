from pylmm.input import plink
from pylmm.lmm import LMM
from os.path import join, split, splitext
import pylmm.pylmmGWAS
from scipy import linalg
import numpy as np
import random


'''
missing phenotype values and missing covariate values, make a copy of data set that we have,
mess with it a lil bit, replace some phenotypes with -9 or NAs
see if they're excluded from the analysis
people who do have missing phenotypes are taken out of the kinship matrix too

assoc tests too,

code for generating data, based on parameters
focus on the missing phenotype masks, try and write/run the gwas by itself,
make sure the K from the readkfile is the same from


make sure if scramble order of individuals, get the same thing back
1) take a dataset that functions, scramble the order of the data, recover the same results
2) add missing values





def gen_missing_phenos(phenolist):
    newphenos = []
    for pheno in phenolist:
        randint = random.randint(0, 6)
        if randint % 5 == 0:
            newphenos.append([pheno[0], pheno[1], str(-9)])
        else:
            newphenos.append(pheno)
    return newphenos
'''



def checkoutput():
    file1_results = []
    with open('out_clean.noX.test', 'r') as file1:
        for line in file1:
            file1_results.append(line)
    file2_results = []
    with open('out_clean.noX_e.test', 'r') as file2:
        for line in file2:
            file2_results.append(line)
    assert(all([file1_results[i] == file2_results[i] for i in range(len(file1_results))]))
    return


if __name__ == "__main__":
    more_args = ['-v',
                 '--kfile', '../data/snps.132k.clean.noX.pylmm.kin',
                 '--bfile', '../data/snps.132k.clean.noX',
                 '--phenofile', '../data/snps.132k.clean.noX_e.fake.phenos',
                 'out_clean.noX_e.test']
    pylmmGWAS.main(more_args)
    checkoutput()
