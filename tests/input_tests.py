from pylmm.input import plink
import os
import sys





def test1():
    return True


def run_all_tests():
    output1 = test1()
    return

if __name__ == "__main__":
    run_all_tests()

    my_plink_object = plink(fbase='../data/snps.132k.clean.noX', kFile='../data/snps.132k.clean.noX.pylmm.kin',
                            phenoFile='../data/snps.132k.clean.noX.fake.phenos', type='b', normGenotype=True,
                            readKFile=False,fastLMM_kinship=False) ## FIX ME