from pylmm.input import plink
import os
import sys





def test1(plink_object):
    ## Ensure the number of individals == the number of lines in the fam file
    counter = 0
    with open('../data/snps.132k.clean.nox.fam', 'r') as opfile:
        for line in opfile:
            counter += 1
    assert counter == len(plink_object.indivs), \
        "getIndivs failure: counted %s lines, got  %i indivs" % (counter, len(plink_object.indivs))
    return

def test2(plink_object):
    ## Ensure the number of individuals == the number of phenotypes read in
    assert len(plink_object.phenos) == len(plink_object.indivs), \
        "Number of phentotypes %s does not equal number of individuals read %r" % \
        (len(plink_object.phenos), len(plink_object.indivs))
    return

def test3(plink_object):
    ## Check that a kinship file has been read
    assert plink_object.kFile is not None, ("Please specify a kinship matrix")
    return


def run_all_tests(plink_object):
    test1(plink_object)
    test2(plink_object)
    test3(plink_object) ## make readKfile true assertion
    return

if __name__ == "__main__":
    my_plink_object = plink(fbase='../data/snps.132k.clean.noX', kFile='../data/snps.132k.clean.noX.pylmm.kin',
                            phenoFile='../data/snps.132k.clean.noX.fake.phenos', type='b', normGenotype=True,
                            readKFile=False,fastLMM_kinship=False) ## FIX ME
    run_all_tests(my_plink_object)