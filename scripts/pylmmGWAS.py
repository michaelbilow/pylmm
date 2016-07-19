#!/usr/bin/python

# pylmm is a python-based linear mixed-model solver with applications to GWAS
# Copyright (C) 2015  Nicholas A. Furlotte (nick.furlotte@gmail.com)

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pdb
import time
import sys
from optparse import OptionParser, OptionGroup
import numpy as np
from scipy import linalg
from pylmm.lmm import LMM
from pylmm import input
from os import listdir
from os.path import join, split, splitext, exists, isfile
import gzip


def printOutHead():
    out.write("\t".join(["SNP_ID", "BETA", "BETA_SD", "F_STAT", "P_VALUE"]) + "\n")


def outputResult(id, beta, betaSD, ts, ps):
    out.write("\t".join([str(x) for x in [id, beta, betaSD, ts, ps]]) + "\n")


def parse_command_line(parser):
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()

    output_fn = args[0]

    if not options.tfile and not options.bfile and not options.emmaFile:
        # if not options.pfile and not options.tfile and not options.bfile:
        parser.error("You must provide at least one PLINK input file base "
                     "(--tfile or --bfile) or an EMMA formatted file (--emmaSNP).")
    if not options.kfile:
        parser.error("Please provide a pre-computed kinship file")

    # READING PLINK input
    if options.verbose:
        sys.stderr.write("Reading SNP input...\n")
    if options.bfile:
        plink_object = input.plink(options.bfile, type='b',
                                   phenoFile=options.phenoFile,
                                   normGenotype=options.normalizeGenotype)
    elif options.tfile:
        plink_object = input.plink(options.tfile, type='t',
                                   phenoFile=options.phenoFile,
                                   normGenotype=options.normalizeGenotype)
    elif options.emmaFile:
        plink_object = input.plink(options.emmaFile,
                                   type='emma',
                                   phenoFile=options.phenoFile,
                                   normGenotype=options.normalizeGenotype)
    else:
        parser.error("You must provide at least one PLINK input file base")
    return options, args, plink_object


def read_phenotypes(parser, plink_object):
    options, args = parser.parse_args()
    if not isfile(options.phenoFile) and not isfile(options.emmaPheno):
        parser.error("Please provide a phenotype file using the --phenofile or --emmaPHENO argument.")

    # Read the emma phenotype file if provided.
    # Format should be rows are phenotypes and columns are individuals.
    if options.emmaPheno:
        with open(options.emmaPheno, 'r') as emma_phenofile:
            P = []
            for line in emma_phenofile:
                v = line.strip().split()
                p = []
                for x in v:
                    try:
                        p.append(float(x))
                    except ValueError:
                        p.append(np.nan)
                P.append(p)
        plink_object.phenos = np.array(P).T
    return plink_object


def read_covariates(parser, plink_object):
    options, args = parser.parse_args()
    # READING Covariate File
    if options.covfile:
        if options.verbose:
            sys.stderr.write("Reading covariate file...\n")
        raw_covariates = plink_object.getCovariates(options.covfile)
        if options.noMean:
            X0 = raw_covariates
        else:
            X0 = np.hstack([np.ones((plink_object.phenos.shape[0], 1)), raw_covariates])
    elif options.emmaCov:
        if options.verbose:
            sys.stderr.write("Reading covariate file...\n")
        raw_covariates = plink_object.getCovariatesEMMA(options.emmaCov)
        if options.noMean:
            X0 = raw_covariates
        else:
            X0 = np.hstack([np.ones((plink_object.phenos.shape[0], 1)), raw_covariates])
    else:
        X0 = np.ones((plink_object.phenos.shape[0], 1))

    if np.isnan(X0).sum():
        parser.error("The covariate file {} contains missing values. "
                     "At this time we are not dealing with this case.  "
                     "Either remove those individuals with missing values "
                     "or replace them in some way.".format(options.covfile))
    return X0


def read_kfile(fn, n_indivs, verbose):
    if verbose:
        sys.stderr.write("Reading kinship...\n")
    begin = time.time()
    # This method seems to be the fastest and works if you already know the size of the matrix
    if splitext(fn) == '.gz':
        with gzip.open(fn, 'r') as input_zip:
            K = np.fromstring(input_zip, sep=' ')  # Assume that space separated
    else:
        K = np.fromfile(open(fn, 'r'), sep=" ")
    K.resize((len(n_indivs), len(n_indivs)))

    end = time.time()
    # Other slower ways
    # K = np.loadtxt(options.kfile)
    # K = np.genfromtxt(options.kfile)
    if verbose:
        sys.stderr.write("Read the {} x {} kinship matrix in {:.3} seconds \n"
                         .format(K.shape[0], K.shape[1], end - begin))
    return K


def read_kinship(parser, plink_object):
    options, args = parser.parse_args()
    K = read_kfile(options.kfile, len(plink_object.indivs), options.verbose)
    if options.kfile2:  ## Currently Deprecated
        K2 = read_kfile(options.kfile2, len(plink_object.indivs), options.verbose)
    else:
        K2 = None
    return K, K2


def setup_LMM(parser, plink_object, Y, X0, K, K2):
    options, args = parser.parse_args()

    # PROCESS the phenotype data -- Remove missing phenotype values
    # Keep will now index into the "full" data to select what we keep (either everything or a subset of non missing data
    Y = plink_object.phenos[:, options.pheno]
    v = np.isnan(Y)
    keep = True - v
    if v.sum():
        if options.verbose:
            sys.stderr.write("Cleaning the phenotype vector by removing %d individuals...\n".format(v.sum()))
        Y = Y[keep]
        X0 = X0[keep, :]
        K = K[keep, :][:, keep]
        if options.kfile2:
            K2 = K2[keep, :][:, keep]
        Kva = []
        Kve = []

    # Only load the decomposition if we did not remove individuals.
    # Otherwise it would not be correct and we would have to compute it again.
    if not v.sum() and options.eigenfile:
        if options.verbose:
            sys.stderr.write("Loading pre-computed eigendecomposition...\n")
        Kva = np.load(options.eigenfile + ".Kva")
        Kve = np.load(options.eigenfile + ".Kve")
    else:
        Kva = []
        Kve = []

    # CREATE LMM object for association
    n = K.shape[0]
    if not options.kfile2:
        lmm_object = LMM(Y, K, Kva, Kve, X0, verbose=options.verbose)
    else:
        lmm_object = LMM_withK2(Y, K, Kva, Kve, X0, verbose=options.verbose, K2=K2) # Not implemented

    if options.verbose:
        sys.stderr.write("Computing fit for null model\n")
    lmm_object.fit()

    if options.verbose and not options.kfile2:
        sys.stderr.write("\t** heritability={:.3}, sigma={:.3}\n".format(lmm_object.optH, lmm_object.optSigma))
    if options.verbose and options.kfile2: sys.stderr.write(
        "\t** heritability=%0.3f, sigma=%0.3f, w=%0.3f\n" % (lmm_object.optH, lmm_object.optSigma, lmm_object.optW))
    return lmm_object

# Buffers for pvalues and t-stats
PS = []
TS = []
count = 0
out = open(output_fn, 'w')
printOutHead()

for snp, id in plink_object:
    count += 1
    if options.verbose and count % 1000 == 0:
        sys.stderr.write("At SNP %d\n" % count)

    x = snp[keep].reshape((n, 1))
    v = np.isnan(x).reshape((-1,))
    # Check SNPs for missing values
    if v.sum():
        keeps = True - v
        xs = x[keeps, :]
        if keeps.sum() <= 1 or xs.var() <= 1e-6:
            PS.append(np.nan)
            TS.append(np.nan)
            outputResult(id, np.nan, np.nan, np.nan, np.nan)
            continue

        # Its ok to center the genotype -  I used options.normalizeGenotype to
        # force the removal of missing genotypes as opposed to replacing them with MAF.
        if not options.normalizeGenotype: xs = (xs - xs.mean()) / np.sqrt(xs.var())
        Ys = Y[keeps]
        X0s = X0[keeps, :]
        Ks = K[keeps, :][:, keeps]
        if options.kfile2: K2s = K2[keeps, :][:, keeps]
        if options.kfile2:
            Ls = LMM_withK2(Ys, Ks, X0=X0s, verbose=options.verbose, K2=K2s)
        else:
            Ls = LMM(Ys, Ks, X0=X0s, verbose=options.verbose)
        if options.refit:
            Ls.fit(X=xs, REML=options.REML)
        else:
            # try:
            Ls.fit(REML=options.REML)
            # except: pdb.set_trace()
        ts, ps, beta, betaVar = Ls.association(xs, REML=options.REML, returnBeta=True)
    else:
        if x.var() == 0:
            PS.append(np.nan)
            TS.append(np.nan)
            outputResult(id, np.nan, np.nan, np.nan, np.nan)
            continue

        if options.refit: L.fit(X=x, REML=options.REML)
        ts, ps, beta, betaVar = L.association(x, REML=options.REML, returnBeta=True)

    outputResult(id, beta, np.sqrt(betaVar).sum(), ts, ps)
    PS.append(ps)
    TS.append(ts)

if __name__ == "__main__":
    pass