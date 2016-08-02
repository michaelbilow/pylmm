from optparse import OptionParser, OptionGroup
import sys
from os.path import exists, isfile

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

def ALBI_parser():
    usage = """usage: %prog [options] --albi_kinship_eigenvalues filename --albi_kinship_eigenvectors filename
                --albi_covariates filename --albi_save_dist_filename filename
    ALBI can be used to estimate the distributions of heritability and their confidence intervals and to save the results.
    There are two cases in which ALBI distributions can be made:

        1) All covariates are eigenvectors:
        python albi.py --albi_kinship_eigenvalues filename
                 [--albi_use_eigenvectors_as_covariates <list of #>]
                 [--albi_precision <# of grid points>]
                 [--albi_distribution_precision <# of grid points>]
                 [--albi_samples <# of random samples>]
                  --albi_save_dist_filename filename

        2) General covariates
        python albi.py --albi_kinship_eigenvalues filename
                  --albi_kinship_eigenvectors filename
                  --albi_covariates filename
                 [--albi_precision <# of grid points>]
                 [--albi_distribution_precision <# of grid points>]
                 [--albi_samples <# of random samples>]
                  --albi_save_dist_filename filename

        To load ALBI estimated distributions and calculate their confidence intervals
        python albi.py --albi_load_dist_filename filename
                 (--albi_estimates_filename filename
                    or
                  --albi_estimate_grid <# of grid points>)
                 [--albi_confidence <required confidence level>]
                 [--albi_output_filename filename]
                """
    parser = OptionParser(usage=usage)
    eigenGroup = OptionGroup(parser, "Covariates are eigenvectors")
    gencovGroup = OptionGroup(parser, "General covariates")
    CIGroup = OptionGroup(parser, "Create CIs")

    eigenGroup.add_option("-albi_k", "--albi_kinship_eigenvalues", dest='albi_eigenvalfile',
                          help='A file containing the eigenvalues of the kinship matrix, one '
                               'eigenvalue per line, in text format. This could be created, for '
                               'example, with GCTAs --pca flag.')
    eigenGroup.add_option('-albi_u', '--albi_use_eigenvectors_as_covariates', dest='albi_eigenvec_list',
                          help='A list detailing which eigenvectors should be used as covariates. For '
                               'example: If you only use an intercept term, and your kinship matrix was '
                               'constructed from a standardized SNP matrix, use -1 (the last eigenvector '
                               'is equal to the constant vector); if in addition, you add the first 3 PCs '
                               'as covariates, use 0,1,2,-1. Default is -1.')
    eigenGroup.add_option('-albi_p', '--albi_precision', dest='albi_precision',
                          help='he number of grid points of the true heritability values, for which the '
                               'estimator distributions are estimated. Effectively, this is the precision '
                               'at which the CIs will be given (e.g., 100 grid points = 0.01 precision). '
                               'Default is 100.')
    eigenGroup.add_option('-albi_d', '--albi_distribution_precision', dest='albi_dist_precision',
                          help='The number of grid points at which each estimator distribution is estimated. '
                               'This controls the accuracy of estimation. Default is 100.')
    eigenGroup.add_option('-albi_n', '--albi_samples', dest='albi_samples',
                          help='Filename at which to save the estimated distributions.')
    eigenGroup.add_option('-albi_s', '--albi_dest_filename', dest='albi_savefile',
                           help='Filename at which to save the estimated distributions.')

    gencovGroup.add_option('-albi_k', '--albi_kinship_eigenvalues', dest='albi_eigenvalfile',
                           help='A file containing the eigenvalues of the kinship matrix, one '
                                'eigenvalue per line, in text format. This could be created, for '
                                'example, with GCTAs --pca flag.')
    gencovGroup.add_option('-albi_v', '--albi_kinship_eigenvectors', dest='albi_eigenvecfile',
                           help='A file containing the eigenvectors of the kinship matrix, one '
                                'eigenvector per column, in text format. This could be created, '
                                'for example, with GCTAs --pca flag.')
    gencovGroup.add_option('-albi_x', '--albi_covariates', dest='albi_covfile',
                           help='A file containing the covariates, one covariate per column, in text format.'
                                ' Remember to include a constant column if you used an intercept term.')
    gencovGroup.add_option('-albi_p', '--albi_precision', dest='albi_precision',
                           help='The number of grid points of the true heritability values, for which the '
                                'estimator distributions are estimated. Effectively, this is the precision '
                                'at which the CIs will be given (e.g., 100 grid points = 0.01 precision). '
                                'Default is 100.')
    gencovGroup.add_option('-albi_d', '--albi_distribution_precision', dest='albi_dist_precision',
                           help='The number of grid points at which each estimator distribution is estimated.'
                                ' This controls the accuracy of estimation. Default is 100.')
    gencovGroup.add_option('-albi_n', '--albi_samples', dest='albi_samples',
                           help='Number of random bootstrap samples to use for estimation. Default is 1000.')
    gencovGroup.add_option('-albi_s', '--albi_dest_filename', dest='albi_savefile',
                           help='Filename at which to save the estimated distributions.')

    CIGroup.add_option('-albi_l', '--albi_load_dist_filename', dest='albi_loadfile',
                       help='Filename from which to load the estimated distributions.')
    CIGroup.add_option('-albi_f', '--albi_estimates_filename', dest='albi_estfile',
                       help='A filename containing a list of heritability estimates '
                            '(one per line) in text format. A CI will be calculated for each one.')
    CIGroup.add_option('-albi_g', '--albi_estimate_grid', dest='estgrid',
                       help='Alternatively, one can ask ALBI to calculate CIs for a grid of heritability '
                            'estimates (e.g., a grid of 100, will calculate CIs for 0, 0.01, ..., 0.99, 1). '
                            'Default is 10.')
    CIGroup.add_option('-albi_c', '--albi_confidence', dest='albi_conf',
                       help='The required confidence level for the CIs. Default is 0.95 (95% CIs).')
    CIGroup.add_option('-albi_o', '--albi_output_filename', dest='albi_output',
                       help='File to which to write the calculated CIs. If not supplied, CIs will be printed.')

    parser.add_option_group(eigenGroup)
    parser.add_option_group(gencovGroup)
    parser.add_option_group(CIGroup)

    options, args = parser.parse_args()
    return options, args

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


    options, args = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()

    # Parser Assetions

    if not options.tfile and not options.bfile and not options.emmaFile:
        # if not options.pfile and not options.tfile and not options.bfile:
        parser.error("You must provide at least one PLINK input file base "
                     "(--tfile or --bfile) or an EMMA formatted file (--emmaSNP).")
    if not options.kfile:
        parser.error("Please provide a pre-computed kinship file")

    if not options.bfile and not options.tfile and not options.emmaFile:
        parser.error("You must provide at least one PLINK input file base")

    # if (not options.phenoFile and not options.emmaPheno) or not isfile(options.phenoFile) and not isfile(options.emmaPheno):
    #     print options.phenoFile
    #     print options.emmaPheno
    #     parser.error("Please provide a phenotype file using the --phenofile or --emmaPHENO argument.")

    return options, args
