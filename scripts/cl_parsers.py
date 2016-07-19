from optparse import OptionParser, OptionGroup


def kinship_parser():
    """
    Sets up a command line parser for pylmmKinship.py
    """
    usage = """usage: %prog [options] --[tfile | bfile] plinkFileBase outfile

    """
    parser = OptionParser(usage=usage)

    basicGroup = OptionGroup(parser, "Basic Options")
    # advancedGroup = OptionGroup(parser, "Advanced Options")

    # basicGroup.add_option("--pfile", dest="pfile",
    #                  help="The base for a PLINK ped file")
    basicGroup.add_option("--tfile", dest="tfile",
                          help="The base for a PLINK tped file")
    basicGroup.add_option("--bfile", dest="bfile",
                          help="The base for a PLINK binary ped file")
    basicGroup.add_option("--emmaSNP", dest="emmaFile", default=None,
                          help="For backwards compatibility with emma, we allow for \"EMMA\" file formats."
                               "  This is just a text file with individuals on the rows and snps on the columns.")
    basicGroup.add_option("--emmaNumSNPs", dest="numSNPs", type="int", default=0,
                          help="When providing the emmaSNP file you need to specify how many snps are in the file")

    basicGroup.add_option("-e", "--efile", dest="saveEig", help="Save eigendecomposition to this file.")
    basicGroup.add_option("-n", default=1000, dest="computeSize", type="int",
                          help="The maximum number of SNPs to read into memory at once (default 1000).  "
                               "This is important when there is a large number of SNPs, "
                               "because memory could be an issue.")

    basicGroup.add_option("-v", "--verbose",
                          action="store_true", dest="verbose", default=False,
                          help="Print extra info")

    parser.add_option_group(basicGroup)
    # parser.add_option_group(advancedGroup)
    return parser


def GWAS_parser():
    usage = """usage: %prog [options] --kfile kinshipFile --[tfile | bfile] plinkFileBase outfileBase

    This program provides basic genome-wide association (GWAS) functionality.
    You provide a phenotype and genotype file as well as a pre-computed (use pylmmKinship.py)
    kinship matrix and the program outputs a result file with information about each SNP,
    including the association p-value.
    The input file are all standard plink formatted with the first two columns specifiying
    the individual and family ID.  For the phenotype file, we accept either NA or -9 to denote missing values.

    Basic usage:

        python pylmmGWAS.py -v --bfile plinkFile --kfile preComputedKinship.kin --phenofile plinkFormattedPhenotypeFile resultFile

    """
    parser = OptionParser(usage=usage)

    basicGroup = OptionGroup(parser, "Basic Options")
    advancedGroup = OptionGroup(parser, "Advanced Options")
    experimentalGroup = OptionGroup(parser, "Experimental Options")

    # basicGroup.add_option("--pfile", dest="pfile",
    #                  help="The base for a PLINK ped file")
    basicGroup.add_option("--tfile", dest="tfile",
                          help="The base for a PLINK tped file")
    basicGroup.add_option("--bfile", dest="bfile",
                          help="The base for a PLINK binary bed file")
    basicGroup.add_option("--phenofile", dest="phenoFile", default=None,
                          help="Without this argument the program will look for a file with "
                               ".pheno that has the plinkFileBase root.  If you want to specify an alternative "
                               "phenotype file, then use this argument.  This file should be in plink format. ")

    # EMMA Options
    basicGroup.add_option("--emmaSNP", dest="emmaFile", default=None,
                          help="For backwards compatibility with emma, we allow for \"EMMA\" file formats.  "
                               "This is just a text file with individuals on the columns and snps on the rows.")
    basicGroup.add_option("--emmaPHENO", dest="emmaPheno", default=None,
                          help="For backwards compatibility with emma, we allow for \"EMMA\" file formats.  "
                               "This is just a text file with each phenotype as one row.")
    basicGroup.add_option("--emmaCOV", dest="emmaCov", default=None,
                          help="For backwards compatibility with emma, we allow for \"EMMA\" file formats.  "
                               "This is just a text file with each covariate as one row.")

    basicGroup.add_option("--kfile", dest="kfile",
                          help="The location of a kinship file.  This is an nxn plain text file and can be "
                               "computed with the pylmmKinship program.")
    basicGroup.add_option("--covfile", dest="covfile",
                          help="The location of a covariate file file.  This is a plink formatted covariate file.")
    basicGroup.add_option("-p", type="int", dest="pheno",
                          help="The phenotype index to be used in association.", default=0)

    advancedGroup.add_option("--removeMissingGenotypes",
                             action="store_false", dest="normalizeGenotype", default=True,
                             help="By default the program replaces missing genotypes with the minor allele "
                                  "frequency.  This option overrides that behavior making the program remove "
                                  "missing individuals.  NOTE: This can increase running time due to the need to "
                                  "recompute the eigendecomposition for each SNP with missing values.")
    advancedGroup.add_option("--refit",
                             action="store_true", dest="refit", default=False,
                             help="Refit the variance components at each SNP (default is to lock in the "
                                  "variance components under the null).")

    advancedGroup.add_option("--REML",
                             action="store_true", dest="REML", default=False,
                             help="Use restricted maximum-likelihood (REML) (default is maximum-likelihood).")
    # advancedGroup.add_option("-e", "--efile", dest="saveEig", help="Save eigendecomposition to this file.")
    advancedGroup.add_option("--eigen", dest="eigenfile",
                             help="The location of the precomputed eigendecomposition for the kinship file.  "
                                  "These can be computed with pylmmKinship.py.")
    advancedGroup.add_option("--noMean", dest="noMean", default=False, action="store_true",
                             help="This option only applies when --covfile is used.  When covfile is provided, "
                                  "the program will automatically add a global mean covariate to the model "
                                  "unless this option is specified.")

    advancedGroup.add_option("-v", "--verbose",
                             action="store_true", dest="verbose", default=False,
                             help="Print extra info")

    # Experimental Group Options
    experimentalGroup.add_option("--kfile2", dest="kfile2",
                                 help="The location of a second kinship file.  This file has the same "
                                      "format as the first kinship.  This might be used if you want to "
                                      "correct for another form of confounding.")

    parser.add_option_group(basicGroup)
    parser.add_option_group(advancedGroup)
    parser.add_option_group(experimentalGroup)
