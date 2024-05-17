#!/usr/bin/env python

from Bio.Seq import Seq
import argparse
import gzip
import itertools
import json
import logging
import os
import re
import sys
import pandas as pd
import pyfastx


OUTS_DIR = 'outs'
FILTERED_MATRIX_DIR_SUFFIX = 'filtered'
BARCODE_FILE_NAME = 'barcodes.tsv.gz'

logger = logging.getLogger(__name__)


def openfile(file_name, mode='rt', **kwargs):
    """open gzip or plain file"""
    if file_name.endswith('.gz'):
        file_obj = gzip.open(file_name, mode=mode, **kwargs)
    else:
        file_obj = open(file_name, mode=mode, **kwargs)
    return file_obj


def get_matrix_file_path(matrix_dir, file_name):
    """
    compatible with non-gzip file
    """
    non_gzip = file_name.strip('.gz')
    file_path_list = [f'{matrix_dir}/{file_name}', f'{matrix_dir}/{non_gzip}']
    for file_path in file_path_list:
        if os.path.exists(file_path):
            return file_path
        
        
def get_matrix_dir_from_match_dir(match_dir):
    """
    Returns:
        matrix_dir: PosixPath object
    """
    matrix_dir = f"{match_dir}/{OUTS_DIR}/{FILTERED_MATRIX_DIR_SUFFIX}"
    if not os.path.exists(matrix_dir):
        raise FileNotFoundError(f'{matrix_dir} not found')
    
    return matrix_dir


def get_barcode_from_matrix_dir(matrix_dir):
    """
    Returns:
        match_barcode: list
        no_match_barcode: int
    """
  
    match_barcode_file = get_matrix_file_path(matrix_dir, BARCODE_FILE_NAME)
    match_barcode, n_match_barcode = read_one_col(match_barcode_file)

    return match_barcode, n_match_barcode


def get_barcode_from_match_dir(match_dir):
    '''
    multi version compatible
    Returns:
        match_barcode: list
        no_match_barcode: int
    '''
    matrix_dir = get_matrix_dir_from_match_dir(match_dir)
    return get_barcode_from_matrix_dir(matrix_dir)


def get_seq_str(seq, sub_pattern):
    """
    join seq slices.

    Args:
        seq: usually R1 read
        sub_pattern: [slice(0,8),slice(16,24)]

    Returns:
        joined intervals seq

    >>> sub_pattern_dict = [slice(0,2)]
    >>> seq = "A" * 2 + "T" * 2
    >>> get_seq_str(seq, sub_pattern_dict)
    'AA'
    """
    return "".join([seq[x] for x in sub_pattern])


def get_seq_list(seq, pattern_dict, abbr):
    """
    >>> pattern_dict = Barcode.parse_pattern("C2L3C2")
    >>> seq = "AAGGGTT"
    >>> Barcode.get_seq_list(seq, pattern_dict, "C")
    ['AA', 'TT']
    """
        
    return [seq[item[0]: item[1]] for item in pattern_dict[abbr]]


def findall_mismatch(seq, n_mismatch=1, bases='ACGTN'):
    """
    choose locations where there's going to be a mismatch using combinations
    and then construct all satisfying lists using product

    Return:
    all mismatch <= n_mismatch set.

    >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
    >>> seq_set = findall_mismatch("ACG")
    >>> seq_set == answer
    True
    """
    seq_set = set()
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for locs in itertools.combinations(range(seq_len), n_mismatch):
        seq_locs = [[base] for base in seq]
        for loc in locs:
            seq_locs[loc] = list(bases)
        for poss in itertools.product(*seq_locs):
            seq_set.add(''.join(poss))
    return seq_set

def get_mismatch_dict(seq_list, n_mismatch=1):
    """
    Return:
    mismatch dict. Key: mismatch seq, value: seq in seq_list

    >>> seq_list = ["AACGTGAT", "AAACATCG"]
    >>> mismatch_dict = get_mismatch_dict(seq_list)
    >>> mismatch_dict["AACGTGAA"] == "AACGTGAT"
    True
    """
    mismatch_dict = {}
    for seq in seq_list:
        seq = seq.strip()
        if seq == '':
            continue
        for mismatch_seq in findall_mismatch(seq, n_mismatch):
            mismatch_dict[mismatch_seq] = seq
    return mismatch_dict


def read_one_col(file):
    """
    Read file with one column. Strip each line.
    Returns col_list, line number
    """
    df = pd.read_csv(file, header=None)
    col1 = list(df.iloc[:, 0])
    col1 = [item.strip() for item in col1]
    num = len(col1)
    return col1, num


def parse_pattern(pattern, allowed="CLUNT"):
    """
    >>> pattern_dict = parse_pattern("C8L16C8L16C8L1U12T18")
    >>> pattern_dict['C']
    [slice(0, 8, None), slice(24, 32, None), slice(48, 56, None)]
    >>> pattern_dict['L']
    [slice(8, 24, None), slice(32, 48, None), slice(56, 57, None)]
    """
    pattern_dict = {}
    p = re.compile(r'([A-Z])(\d+)')
    tmp = p.findall(pattern)
    if not tmp:
        sys.exit(f'Invalid pattern: {pattern}')
    start = 0
    for x, length in tmp:
        if x not in allowed:
            sys.exit(f'Invalid pattern: {pattern}')
        if x not in pattern_dict:
            pattern_dict[x] = []
        end = start + int(length)
        pattern_dict[x].append(slice(start,end))
        start = end
    return pattern_dict


def get_raw_mismatch(files: list, n_mismatch: int):
    """
    Args:
        files: whitelist file paths
        n_mismatch: allowed number of mismatch bases
    Returns:
        raw_list
        mismatch_list
    """    
    raw_list, mismatch_list = [], []
    for f in files:
        barcodes, _ = read_one_col(f)
        raw_list.append(set(barcodes))
        barcode_mismatch_dict = get_mismatch_dict(barcodes, n_mismatch)
        mismatch_list.append(barcode_mismatch_dict)

    return raw_list, mismatch_list


def check_seq_mismatch(seq_list, raw_list, mismatch_list):
    '''
    Returns
        valid: True if seq in mismatch_list
        corrected: True if seq in mismatch_list but not in raw_list
        res: joined seq

    >>> seq_list = ['ATA', 'AAT', 'ATA']
    >>> correct_set_list = [{'AAA'},{'AAA'},{'AAA'}]
    >>> mismatch_dict_list = [get_mismatch_dict(['AAA'])] * 3

    >>> check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
    (True, True, 'AAA_AAA_AAA')

    >>> seq_list = ['AAA', 'AAA', 'AAA']
    >>> check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
    (True, False, 'AAA_AAA_AAA')
    '''
    valid = True
    corrected = False
    res = []
    for index, seq in enumerate(seq_list):
        if seq not in raw_list[index]:
            if seq not in mismatch_list[index]:
                valid = False
                res = []
            else:
                corrected = True
                res.append(mismatch_list[index][seq])
        else:
            res.append(seq)

    return valid, corrected, '_'.join(res)

def get_protocol_dict(assets_dir):
    """
    Return:
    protocol_dict. Key: protocol name, value: protocol dict

    >>> protocol_dict = get_protocol_dict("./assets/")
    >>> protocol_dict["GEXSCOPE-MicroBead"]["pattern_dict"]
    {'C': [slice(0, 12, None)], 'U': [slice(12, 20, None)]}
    """
    json_file = os.path.join(assets_dir, "protocols.json")
    protocol_dict = json.load(open(json_file))
    whitelist_dir = os.path.join(assets_dir, "whitelist")
    # add folder prefix
    for protocol in protocol_dict:
        cur = protocol_dict[protocol]
        bc = cur.get("bc", [])
        linker = cur.get("linker", [])
        if bc:
            cur["bc"] = [os.path.join(whitelist_dir, protocol, x) for x in bc]
        if linker:
            cur["linker"] = [os.path.join(whitelist_dir, protocol, x) for x in linker]
        cur["pattern_dict"] = parse_pattern(cur["pattern"])
    return protocol_dict


def reverse_complement(seq):
    """Reverse complementary sequence

    :param original seq
    :return Reverse complementary sequence
    """
    return str(Seq(seq).reverse_complement())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--fq1', required=True)
    parser.add_argument('--fq2', required=True)
    parser.add_argument('--assets_dir', required=True)
    parser.add_argument('--protocol', required=True)
    parser.add_argument('--match_dir', required=True)
    args = parser.parse_args()
    
    # metrics 
    raw_reads = 0
    valid_reads = 0
    valid_matched_reads = 0
    valid_matched_barcodes = set()
    match_barcodes = set(get_barcode_from_match_dir(args.match_dir)[0]) # barcode set of flv_rna.

    # protocol
    protocol_dict = get_protocol_dict(args.assets_dir)
    protocol = protocol_dict[args.protocol]
    pattern_dict = protocol["pattern_dict"]
    raw_list,  mismatch_list = get_raw_mismatch(protocol["bc"], 1)
        
    # out_fq
    out_fq_fn = {x: f"{args.sample}_R{x}.fq.gz" for x in [1,2]}
    outdict = {k:openfile(v,'wt') for k,v in out_fq_fn.items()}

    fq1_list = args.fq1.split(',')
    fq2_list = args.fq2.split(',')
    raw_reads = 0
    for fq1,fq2 in zip(fq1_list, fq2_list):
        fq1 = pyfastx.Fastx(fq1)
        fq2 = pyfastx.Fastx(fq2)

        for (name1, seq1, qual1), (name2,seq2,qual2) in zip(fq1, fq2):
            raw_reads += 1
            seq_list = get_seq_list(seq1, pattern_dict, 'C')
            # flv
            seq_list = [reverse_complement(seq) for seq in seq_list[::-1]]
            valid, corrected, corrected_seq = check_seq_mismatch(seq_list, raw_list, mismatch_list)
            
            if valid:
                valid_reads += 1
                
            umi = get_seq_str(seq1, pattern_dict['U'])  
            if not umi:
                continue              
            if corrected_seq in match_barcodes:
                valid_matched_reads += 1
                valid_matched_barcodes.add(corrected_seq)
                qual1 = 'F' * len(corrected_seq + umi)
                outdict[2].write(f'@{corrected_seq}:{umi}:{raw_reads}\n{seq2}\n+\n{qual2}\n')
                outdict[1].write(f'@{corrected_seq}:{umi}:{raw_reads}\n{corrected_seq}{umi}\n+\n{qual1}\n')

    data_dict = {
        "sample": args.sample,
        "Raw Reads": raw_reads,
        "Valid Reads": valid_reads,
        "Valid Match Reads": valid_matched_reads,
        "Matched Barcodes": len(valid_matched_barcodes),
    }
    
    stats_file = args.sample + ".protocol_stats.json"
    with open(stats_file, "w") as f:
        json.dump(data_dict, f)