from pylmm.input import plink
from pylmm.lmm import LMM
from os.path import join, split, splitext
from scipy import linalg
import numpy as np
import sys

def lmm_test1(lmm_object, plink_object):
    ## Test that the calculated eigen value is the same as that from the kinship matrix
    kin = []
    with open(plink_object.kFile, 'r') as opfile:
        for line in opfile:
            kin.append([float(x) for x in line.strip().split()])
    kin = np.array(kin)
    for i in range(kin.shape[0]):
        for j in range(kin.shape[1]):
            assert lmm_object.K[i][j] == kin[i][j], "Kinship file not read properly"
    assert np.isclose(lmm_object.Kva, linalg.eigh(kin)[0], atol=1e-05).all() and np.isclose(lmm_object.Kve, linalg.eigh(kin)[1], atol=1e-05).all()
    return





def run_all_lmm_tests(lmm_object, plink_object):
    lmm_test1(lmm_object, plink_object)
    return

if __name__ == "__main__":
    input_folder = '../data/'
    mouse_data = True
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

    lmm_object = LMM(my_plink_object.phenos, my_plink_object.K, verbose=True)
    run_all_lmm_tests(lmm_object, my_plink_object)