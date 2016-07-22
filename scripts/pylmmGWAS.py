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
from cl_parsers import GWAS_parser
import numpy as np
from pylmm.lmm import LMM
from pylmm import input
from os import listdir
from os.path import join, split, splitext, exists, isfile
import gzip


def read_snp_input(options):
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
        raise ValueError("No SNP File!")
    return plink_object


def read_emma_phenos(options, plink_object):
    # Read the emma phenotype file if provided.
    # Format should be rows are phenotypes and columns are individuals.
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


def read_covariates(options, plink_object):
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
        raise ValueError('The covariate file {} contains missing values. '
                         'At this time we are not dealing with this case.  \n'
                         'Either remove those individuals with missing values '
                         'or replace them in some way.'.format(options.covfile))
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
        K = np.fromfile(open(fn, 'r'), sep=' ')
    K.resize(n_indivs, n_indivs)

    end = time.time()
    # Other slower ways
    # K = np.loadtxt(options.kfile)
    # K = np.genfromtxt(options.kfile)
    if verbose:
        sys.stderr.write("Read the {} x {} kinship matrix in {:.3} seconds\n"
                         .format(K.shape[0], K.shape[1], end - begin))
    return K


def read_kinship(options, plink_object):
    K = read_kfile(options.kfile, len(plink_object.indivs), options.verbose)
    if options.kfile2:  ## Currently Deprecated
        K2 = read_kfile(options.kfile2, len(plink_object.indivs), options.verbose)
    else:
        K2 = None
    return K, K2


def load_phenos_and_identify_missing(options, plink_object):
    # PROCESS the phenotype data -- Remove missing phenotype values
    # Keep will now index into the "full" data to select what we keep (either everything or a subset of non missing data
    Y = plink_object.phenos[:, options.pheno]
    missing_pheno_mask = np.isnan(Y)
    return missing_pheno_mask, Y


def clean_missing_phenos(options, missing_pheno_mask, Y, X0, K, K2):
    keep = ~missing_pheno_mask
    if missing_pheno_mask.sum():
        if options.verbose:
            sys.stderr.write("Cleaning the phenotype vector by removing {} individuals...\n"
                             .format(missing_pheno_mask.sum()))
        Y = Y[keep]
        X0 = X0[keep, :]
        K = K[keep, :][:, keep]
        if options.kfile2:
            K2 = K2[keep, :][:, keep]
        Kva, Kve = [], []
    else:
        if options.eigenfile:
            if options.verbose:
                sys.stderr.write("Loading pre-computed eigendecomposition...\n")
            Kva = np.load(options.eigenfile + ".Kva")
            Kve = np.load(options.eigenfile + ".Kve")
        else:
            Kva, Kve = [], []

    return Y, X0, K, K2, Kva, Kve


def setup_lmm_object(options, Y, K, Kva, Kve, X0):
    # CREATE LMM object for association
    if not options.kfile2:
        lmm_object = LMM(Y, K, Kva, Kve, X0, verbose=options.verbose)
    else:  # Not implemented
        # lmm_object = LMM_withK2(Y, K, Kva, Kve, X0, verbose=options.verbose, K2=K2)
        raise ValueError("2 Kinship Matrices Not Implemented; use GCTA to compute variance components")

    if options.verbose:
        sys.stderr.write("Computing fit for null model\n")
    lmm_object.fit()

    if options.verbose and not options.kfile2:
        sys.stderr.write("\t** heritability={:.3}, sigma={:.3}\n".format(lmm_object.optH, lmm_object.optSigma))
    if options.verbose and options.kfile2: sys.stderr.wrte(
        "\t** heritability={:.3}, sigma={:.3}, w={:.3}\n".format(lmm_object.optH, lmm_object.optSigma, lmm_object.optW))
    return lmm_object


def main(fake_args=None):
    print sys.argv
    if fake_args:
        sys.argv += fake_args
    options, args = GWAS_parser()
    plink_object = read_snp_input(options)
    if options.emmaPheno:
        read_emma_phenos(options, plink_object)

    X0 = read_covariates(options, plink_object)
    K, K2 = read_kinship(options, plink_object)

    missing_pheno_mask, Y = load_phenos_and_identify_missing(options, plink_object)
    Y, X0, K, K2, Kva, Kve = clean_missing_phenos(options, missing_pheno_mask, Y, X0, K, K2)

    lmm_object = setup_lmm_object(options, Y, K, Kva, Kve, X0)
    
    output_fn = args[0]
    run_association_tests(options, plink_object, lmm_object, missing_pheno_mask, output_fn,
                          Y, X0, K, K2)
    return


def outputResult(output_file, id, beta, betaSD, ts, ps):
    output_file.write("\t".join([str(x) for x in [id, beta, betaSD, ts, ps]]) + "\n")
    return


def run_association_tests(options, plink_object, lmm_object, missing_pheno_mask, output_fn,
                          Y, X0, K, K2):
    count = 0
    keep = ~missing_pheno_mask
    with open(output_fn, 'w') as output_file:
        output_file.write("\t".join(["SNP_ID", "BETA", "BETA_SD", "F_STAT", "P_VALUE"]) + "\n")

        for snp, snp_id in plink_object:
            count += 1
            if options.verbose and count % 1000 == 0:
                sys.stderr.write("At SNP {}\n".format(count))

            snp_vals = snp[keep].reshape((keep.sum(), 1))
            missing_snp_mask = np.isnan(snp_vals).reshape((-1,))
            # Check SNPs for missing values
            if missing_snp_mask.sum():
                keep_snp_mask = ~missing_snp_mask
                tmp_snp_vals = snp_vals[keep_snp_mask, :]
                if keep_snp_mask.sum() <= 1 or tmp_snp_vals.var() <= 1e-6:
                    outputResult(output_file, snp_id, np.nan, np.nan, np.nan, np.nan)
                    continue

                # Its ok to center the genotype -  I used options.normalizeGenotype to
                # force the removal of missing genotypes as opposed to replacing them with MAF.
                if not options.normalizeGenotype:
                    tmp_snp_vals = (tmp_snp_vals - tmp_snp_vals.mean()) / np.sqrt(tmp_snp_vals.var())

                tmp_Y = Y[keep_snp_mask]
                tmp_X0 = X0[keep_snp_mask, :]
                tmp_K = K[keep_snp_mask, :][:, keep_snp_mask]
                if options.kfile2:
                    K2s = K2[keep_snp_mask, :][:, keep_snp_mask]

                if options.kfile2:
                    tmp_lmm_object = LMM_withK2(tmp_Y, tmp_K, X0=tmp_X0, verbose=options.verbose, K2=K2s)
                else:
                    tmp_lmm_object = LMM(tmp_Y, tmp_K, X0=tmp_X0, verbose=options.verbose)

                tmp_lmm_object.fit(X=tmp_snp_vals, REML=options.REML)
                ts, ps, beta, betaVar = tmp_lmm_object.association(tmp_snp_vals, REML=options.REML, returnBeta=True)
            else:
                if snp_vals.var() == 0:
                    outputResult(snp_id, np.nan, np.nan, np.nan, np.nan)
                    continue

                if options.refit:
                    lmm_object.fit(X=snp_vals, REML=options.REML)
                ts, ps, beta, betaVar = lmm_object.association(snp_vals, REML=options.REML, returnBeta=True)

            outputResult(output_file, snp_id, beta, np.sqrt(betaVar).sum(), ts, ps)

if __name__ == "__main__":
    more_args = ['-v',
                 '--kfile', '../data/snps.132k.clean.noX.pylmm.kin',
                 '--bfile', '../data/snps.132k.clean.noX',
                 '--phenofile', '../data/snps.132k.clean.noX.fake.phenos',
                 'out.test']
    main(more_args)
