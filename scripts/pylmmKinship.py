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

import sys
import pdb
import sys
import os
import numpy as np
from scipy import linalg
from pylmm.lmm import calculateKinship
from pylmm import input
from cl_parsers import kinship_cl_parser


def parse_command_line(parser):
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        sys.exit()

    output_fn = args[0]

    if not options.tfile and not options.bfile and not options.emmaFile:
        parser.error(
            "You must provide at least one PLINK input file base (--tfile or --bfile) or an emma formatted file (--emmaSNP).")

    if options.verbose:
        sys.stderr.write("Reading PLINK input...\n")
    if options.bfile:
        plink_object = input.plink(options.bfile, type='b')
    elif options.tfile:
        plink_object = input.plink(options.tfile, type='t')
    # elif options.pfile: IN = input.plink(options.pfile,type='p')
    elif options.emmaFile:
        if not options.numSNPs:
            parser.error("You must provide the number of SNPs when specifying an emma formatted file.")
        plink_object = input.plink(options.emmaFile, type='emma')
    else:
        parser.error(
            "You must provide at least one PLINK input file base (--tfile or --bfile) or an emma formatted file (--emmaSNP).")

    return plink_object, output_fn, options


def make_kinship(plink_object, options):
    n_indivs = len(plink_object.indivs)
    snps_per_pass = options.computeSize
    W = np.ones((n_indivs, snps_per_pass)) * np.nan

    plink_object.getSNPIterator()
    # Annoying hack to get around the fact that it is expensive to determine the number of SNPs in an emma file
    if options.emmaFile:
        plink_object.numSNPs = options.numSNPs
    i = 0
    K = np.zeros((n_indivs))
    while i < plink_object.numSNPs:
        j = 0
        while j < options.computeSize and i < plink_object.numSNPs:
            snp, _id = plink_object.next()
            if snp.var() == 0:
                i += 1
                continue
            W[:, j] = snp

            i += 1
            j += 1
        if j < options.computeSize:
            W = W[:, range(0, j)]

        if options.verbose:
            sys.stderr.write("Processing first %d SNPs\n" % i)
        if K is None:
            try:
                K = linalg.fblas.dgemm(alpha=1., a=W.T, b=W.T, trans_a=True, trans_b=False)  # calculateKinship(W) * j
            except AttributeError:
                K = np.dot(W, W.T)
        else:
            try:
                K_j = linalg.fblas.dgemm(alpha=1., a=W.T, b=W.T, trans_a=True, trans_b=False)  # calculateKinship(W) * j
            except AttributeError:
                K_j = np.dot(W, W.T)
            K = K + K_j

    K *= 1.0/plink_object.numSNPs
    return K

def write_kinship(K, output_fn, options):
    if options.verbose: sys.stderr.write("Saving Kinship file to %s\n" % outFile)
    np.savetxt(output_fn, K)

    if options.saveEig:
        if options.verbose:
            sys.stderr.write("Obtaining Eigendecomposition\n")
        Kva, Kve = linalg.eigh(K)
        if options.verbose:
            sys.stderr.write("Saving eigendecomposition to {}.[kva | kve]\n".format(output_fn))
        np.savetxt(output_fn + ".kva", Kva)
        np.savetxt(output_fn + ".kve", Kve)

if __name__ == "__main__":
    option_parser = kinship_cl_parser()
    plink_obj, output_fn, parsed_options = parse_command_line(parser=option_parser)
    make_kinship(plink_object=plink_obj, options=parsed_options)