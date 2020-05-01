#! /projects/weinstocklab/comp/local/anaconda2/bin/python

from Bio import SeqIO
import os
import sys
import re
import pickle
import gzip
import logging as L
#import pandas as pd
import numpy as np
#import rpy2.robjects as robjects
#from rpy2.robjects import r as R


def snip(filename, extension=None, alt_extension=None,
         strip_path=False):
    '''return prefix of filename.
    If *extension* is given, make sure that filename has the extension.
    If *strip_path* is set to true, the path is stripped from the file name.
    '''
    if extension:
        if filename.endswith(extension):
            root = filename[:-len(extension)]
        elif alt_extension and filename.endswith(alt_extension):
            root = filename[:-len(alt_extension)]
        else:
            raise ValueError("'%s' expected to end in '%s'" %
                             (filename, extension))
    else:
        root, ext = os.path.splitext(filename)

    if strip_path:
        snipped = os.path.basename(root)
    else:
        snipped = root

    return snipped



def summarizeMapping(infile,ref_db, outfile):
    '''Create a summary dataframe out of the cross match results'''
    #infile, ref_db = infiles


    def _parseCrossMatchDiscrepList(infile, max_len, outfile):
        '''Takes a discrepency list output by cross match (infile), the
        maximum length of the sequences in the reference database to which
        reads were mapped, and the name of the outfile.
        Writes a summary of the mismatches in each aligned read to outfile,
        and returns a nested dictionary of mistmatches that can be converted
        into a multi-indexed dataframe.
        '''

        # set up out dict
        out_dict = {}
        outf = open(outfile, 'w')
        outf.write('Score\tSNP\tINS\tDEL\tQuerySeqID\tQueryStart\tQueryEnd\tQueryTail'
                   '\tReverseComp\tTargetSeqID\tTargetStart\tTargetEnd\tTargetTail\n')


        for line in open(infile):
            if line.startswith('ALIGNMENT'):
                line = line.split()
                length = len(line)
                # assign alginment values to separate IDs
                if length == 13:
                    Score, SNP, INS, DEL, QuerySeqID, QueryStart, QueryEnd, \
                        QueryTail, TargetSeqID, TargetStart, TargetEnd, TargetTail = line[1:]
                    ReverseComp = 'F'
                elif length == 14 and line[-1:] == '*':
                    Score, SNP, INS, DEL, QuerySeqID, QueryStart, QueryEnd, \
                        QueryTail, TargetSeqID, TargetStart, TargetEnd, TargetTail = line[1:-1]
                    ReverseComp = 'F'
                elif len(line) == 14:
                    Score, SNP, INS, DEL, QuerySeqID, QueryStart, QueryEnd, \
                        QueryTail, ReverseComp, TargetSeqID, TargetTail, TargetEnd, TargetStart = line[1:]
                elif len(line) == 15 and line[-1:] == '*':
                    Score, SNP, INS, DEL, QuerySeqID, QueryStart, QueryEnd, \
                        QueryTail, ReverseComp, TargetSeqID, TargetTail, TargetEnd, TargetStart = line[1:-1]
                else:
                    raise ValueError('Expecting either a 13 or 14 field line'
                                     ' output from cross_match ALIGNMENT')

                TargetTail = re.sub('\(|\)', '', TargetTail)
                QueryTail = re.sub('\(|\)', '', QueryTail)

                outf.write('\t'.join([Score, SNP, INS, DEL, QuerySeqID, QueryStart,
                                      QueryEnd, QueryTail, ReverseComp, TargetSeqID,
                                      TargetStart, TargetEnd, TargetTail]) + '\n')

                # sanity check
                if os.path.basename(outfile).endswith('_cross_match_summary.txt'):
                    pass
                else:
                    assert (TargetSeqID, QuerySeqID) not in out_dict.keys(), \
                        ("There are multiple matches for the same sequence to the same"
                         " target in cross_match output: %s %s " % (TargetSeqID, QuerySeqID))

                # reset the SNP/INDEL record and set dict entry
                errors = ['']*max_len
                out_dict[(TargetSeqID, QuerySeqID)] = errors

            elif line.startswith('DISCREPANCY'):
                # $2 == type of discrepency, $5 == discrepency location (1-based)
                # WARNING: All indels are marked as I/D... the number of bases is
                # ignored
                disc_loc = int(line.split()[4]) - 1 # table is 0-based
                expr = re.compile('-[0-9]+')
                disc_type = line.split()[1]
                disc_type = re.sub(expr, '', disc_type)


                # update the SNP/INDEL record
                out_dict[(TargetSeqID, QuerySeqID)][disc_loc] = disc_type

            else:
                continue

        return out_dict

    # WARNING!
    if os.path.basename(outfile).endswith('_cross_match_summary.txt'):
        L.warn('There is no sanity check for duplicate matches in'
               ' default cross_match run')
    #ref_db="/projects/weinstock-lab/students/reference_db/Mock36/mock36reference.fa"
    # get the maximum length of the sequences in the reference database
    # (informs output table size)
    max_len = 0
    for fasta in SeqIO.parse(open(ref_db), 'fasta'):
        if len(str(fasta.seq)) > max_len:
            max_len = len(str(fasta.seq))
    print (max_len)

    # parse cross match output
    out_dict = _parseCrossMatchDiscrepList(infile, max_len, outfile)

    # pickle and dump out_dict
    out_pickle = snip(outfile, '.txt') + '.p'
    pickle.dump(out_dict, open(out_pickle, 'wb'))


summarizeMapping(sys.argv[1],sys.argv[2],sys.argv[3])