from collections import defaultdict
import pandas as pd
import numpy as np
import pyfastx
import copy
import json
import argparse


MAX_CELL = 2 * 10**5


"""
## Features
- CDR3 filtering: contain stop condon, length <=5, etc..

- If barcode A's two chains CDR3s are identical to another barcode B, and A's chain abundance is significantly lower than B's, filter A.

- If `--target_cell_barcode` is provided, the UMI counts of all contigs originated from target cells are multiplied by a weight(default: 6.0) to better distinguish signal from background noise. `--target_cell_barcode` comes from the cell type annotation results of the RNA library.

- Cell-calling is similar to the rna cell-calling algorithm.

## Output
- `clonetypes.tsv` High-level description for each clonotype.

- `{sample}_all_contig.csv` High-level and detailed annotation for each contig.

- `{sample}_all_contig.fasta` All assembled contig sequences.

- `{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of {sample}_all_contig.csv.
    
- `{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
"""


class Auto():
    """
    threshold = top {percentile}% cell count / coef
    count is usually UMI count.
    >>> array = [50] * 100 + [30] * 100 + [10] * 100 + [4] * 100
    >>> Auto(array, coef=10).run()
    5
    >>> Auto(array, percentile=70, coef=3).run()
    10
    >>> Auto(array, percentile=50, coef=10, expected_cell_num=100).run()
    5
    >>> Auto([1, 2, 20, 30, 40], expected_cell_num=4, percentile=50, coef=10).run()
    2
    """
    def __init__(self, array, percentile=99, coef=3, expected_cell_num=None, **kwargs):
        self.array = [x for x in array if x > 0 ]
        self.percentile = percentile
        self.coef = int(coef)
        self.expected_cell_num = expected_cell_num
        self.kwargs = kwargs
    
    def run(self):
        array = self.array
        if not array:
            return 1

        if not self.expected_cell_num:
            expected_cell_num = len(array)
        else:
            expected_cell_num = self.expected_cell_num
            if expected_cell_num > len(array):
                print('Warning: expected_cell_num > len(array)')
                expected_cell_num = len(array)
                      
        sorted_counts = sorted(array, reverse=True)
        count_cell_percentile = np.percentile(sorted_counts[:expected_cell_num], self.percentile)
        threshold = int(count_cell_percentile / self.coef)

        return threshold


def target_cell_calling(df_UMI_sum, expected_target_cell_num=3000, target_barcodes=None, weight=6, coef=5, 
    percentile=85, umi_col='umis'):
    """
    Args:
        df_UMI_sum: A dataframe with columns highest umi's contig and UMI.
    
    Returns:
        target_contigs_id: list
    >>> df_UMI_sum = pd.DataFrame({"contig_id": ["A", "B", "C", "D", "E"], "UMI": [1, 2, 1, 30, 40]})
    >>> target_contigs_id = target_cell_calling(df_UMI_sum, expected_target_cell_num=5, percentile=80, coef=5, target_barcodes=["A", "C"])
    >>> target_contigs_id == {'A_1', 'C_1', 'D_1', 'E_1'}
    True
    """
    if target_barcodes != None:
        target_barcodes = {i for i in target_barcodes}
    umi_threshold = Auto(list(df_UMI_sum[umi_col]), expected_cell_num=expected_target_cell_num, coef=coef, percentile=percentile).run()

    # avoid change the original dataframe
    df_temp = df_UMI_sum.copy()
    if target_barcodes:
        df_temp[umi_col] = df_temp.apply(
            lambda row:  row[umi_col] * weight if row['barcode'] in target_barcodes else row[umi_col], axis=1)
             
    target_contigs = set(df_temp.loc[df_temp[umi_col] >= umi_threshold].contig_id)

    return target_contigs


def parse_contig_file(sample, barcode_report, annot_fa):
    """
    Generate all_contig_annotation file.
    Generate all_contig_fasta file.
    """

    df = pd.read_csv(barcode_report)

    df['productive'] = df['full_length']
    contig_set = set(df.contig_id)

    # generate all contig fasta file
    # add length of each contig. 
    len_dict = dict()
    all_fa = open(f'{sample}_all_contig.fasta' , 'w')
    fa = pyfastx.Fastx(annot_fa)
    for read in fa:
        len_dict[read.name] = read.comment.split(' ')[0]
        if read.name in contig_set:
            sequence = read.sequence
            all_fa.write('>' + read.name + '\n' + sequence + '\n')    
    all_fa.close()
    df['length'] = df['contig_id'].apply(len_dict.get)

    return df

 
def cell_calling(df, seqtype, trust_report, expected_target_cell_num, target_barcodes, target_weight, coef):
    """
    Common filtering based on CDR3:
    Filter nonfunctional CDR3(shown 'out_of_frame' in cdr3 report), or CDR3 sequences containing "N" in the nucleotide sequence.
    Keep CDR3aa start with C.
    Keep CDR3aa length >= 5.
    Keep no stop codon in CDR3aa.
    Filter low abundance contigs based on a umi cut-off. 

    Target cell barcodes filtering(option, --target_cell_barcode needed):
    Filter low abundance contigs based on a umi cut-off.
    The umi counts of all contigs originated from B cells are multiplied by a weight to 
    better distinguish signal from background noise.
    """
    df.sort_values(by='umis', ascending=False, inplace=True)
    if seqtype == 'BCR':
        df_chain_heavy = df[df['chain']=='IGH']
        df_chain_light = df[(df['chain']=='IGK') | (df['chain']=='IGL')]
    else:
        df_chain_heavy = df[df['chain'] == 'TRA']
        df_chain_light = df[df['chain'] == 'TRB']
    df_chain_heavy = df_chain_heavy.drop_duplicates(['barcode'])
    df_chain_light = df_chain_light.drop_duplicates(['barcode'])
    df_for_clono = pd.concat([df_chain_heavy, df_chain_light], ignore_index=True)
        
    # Common filtering
    trust_report = pd.read_csv(trust_report, sep='\t')
    correct_cdr3 = set(df_for_clono.cdr3).intersection(set(trust_report.CDR3aa))
    correct_cdr3 = [i for i in correct_cdr3 if i.startswith('C')]
    correct_cdr3 = [i for i in correct_cdr3 if len(i)>=5]
    correct_cdr3 = [i for i in correct_cdr3 if 'UAG' or 'UAA' or 'UGA' not in i]
    df_for_clono = df_for_clono[df_for_clono['cdr3'].isin(correct_cdr3)]
        
    # Filter low abundance contigs based on a umi cut-off
    if seqtype == 'BCR':
        df_chain_heavy = df_for_clono[df_for_clono['chain']=='IGH']
        df_chain_light = df_for_clono[(df_for_clono['chain']=='IGK') | (df_for_clono['chain']=='IGL')]
    else:
        df_chain_heavy = df_for_clono[df_for_clono['chain'] == 'TRA']
        df_chain_light = df_for_clono[df_for_clono['chain'] == 'TRB']

    filtered_congtigs_id = set()
    for _df in [df_chain_heavy, df_chain_light]:
        target_contigs = target_cell_calling(
        _df, 
        expected_target_cell_num=expected_target_cell_num, 
        target_barcodes=target_barcodes,
        weight = target_weight,
        coef = coef
        )
        filtered_congtigs_id = filtered_congtigs_id | target_contigs       
        
    df_for_clono = df_for_clono[df_for_clono.contig_id.isin(filtered_congtigs_id)]
        
    df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
    cell_barcodes, filtered_contig = set(df_for_clono_pro['barcode']), set(df_for_clono_pro['contig_id'])

    return df_for_clono, cell_barcodes, filtered_contig 


def filter_fasta(sample, cell_barcodes):
    """Filter all contig fasta file by barcodes which are identified to be cell.

    :param cell_barcodes: all barcodes identified to be cell.
    """
    all_contig_fasta = f'{sample}_all_contig.fasta'
    filter_contig_fasta = f'{sample}_filtered_contig.fasta'

    filter_contig_fasta = open(filter_contig_fasta,'w')
    fa = pyfastx.Fastx(all_contig_fasta)
    for read in fa:
        name = read.name
        cb = '_'.join(name.split('_')[:-1]) # remove the suffix num in "bc1_bc2_bc3_num"
        sequence = read.sequence
        if cb in cell_barcodes:
            filter_contig_fasta.write('>' + name + '\n' + sequence + '\n')

    filter_contig_fasta.close()


def parse_clonotypes(sample, df, df_for_clono, cell_barcodes, filtered_contig):
    """Parse clonotypes from CDR3 and manually add clonotype id for each contig.

    :param df: original contig file.
    :param df_for_clono: contig info after filter.
    :param cell_barcodes: all barcodes identified to be cell.
    :return df_filter_contig: filtered contigs by cell barcodes.
    """
    df_for_clono_pro = df_for_clono[df_for_clono['productive']==True].copy()
    df_for_clono_pro['chain_cdr3aa'] = df_for_clono_pro.loc[:, ['chain', 'cdr3']].apply(':'.join, axis=1)
    df_for_clono_pro['chain_cdr3nt'] = df_for_clono_pro.loc[:,['chain', 'cdr3_nt']].apply(':'.join, axis=1)

    cbs = set(df_for_clono_pro['barcode'])
    clonotypes = open('clonotypes.csv', 'w')
    clonotypes.write('barcode\tcdr3s_aa\tcdr3s_nt\n')
    for cb in cbs:
        temp = df_for_clono_pro[df_for_clono_pro['barcode']==cb]
        temp = temp.sort_values(by='chain', ascending=True)
        aa_chain = ';'.join(list(temp['chain_cdr3aa']))
        nt_chain = ';'.join(list(temp['chain_cdr3nt']))
        clonotypes.write(f'{cb}\t{aa_chain}\t{nt_chain}\n')
    clonotypes.close() 

    df_clonotypes = pd.read_csv('clonotypes.csv', sep='\t', index_col=None)
    contig_with_clonotype = copy.deepcopy(df_clonotypes)
    df_dict = df_clonotypes[["cdr3s_nt", "cdr3s_aa"]].set_index("cdr3s_nt").to_dict(orient='dict')['cdr3s_aa']
    df_clonotypes = df_clonotypes.groupby('cdr3s_nt', as_index=False).agg({'barcode': 'count'})
    df_clonotypes.rename(columns={'barcode': 'frequency'}, inplace=True)
    sum_f = df_clonotypes['frequency'].sum()

    df_clonotypes['proportion'] = df_clonotypes['frequency'].apply(lambda x: x/sum_f)
    df_clonotypes.sort_values(by='frequency', ascending=False, inplace=True)
    df_clonotypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_clonotypes.shape[0]+1)]
    df_clonotypes['cdr3s_aa'] = df_clonotypes['cdr3s_nt'].apply(lambda x:df_dict[x])
    df_clonotypes = df_clonotypes.reindex(columns=['clonotype_id', 'frequency', 'proportion', 'cdr3s_aa', 'cdr3s_nt'])
    df_clonotypes.to_csv('clonotypes.csv', sep=',', index=False) 
    used_for_merge = df_clonotypes[['cdr3s_nt','clonotype_id']]

    df_merge = pd.merge(used_for_merge, contig_with_clonotype, on='cdr3s_nt', how='outer')
    df_merge = df_merge[['barcode', 'clonotype_id']]
    df_all_contig = pd.merge(df_merge, df, on='barcode',how='outer')
    df_all_contig.fillna('',inplace = True)
    df_all_contig = df_all_contig[['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'clonotype_id']]
    df_filter_contig = df_all_contig[df_all_contig['barcode'].isin(cell_barcodes)]
    for _df in [df_all_contig, df_filter_contig]:
        _df.loc[~_df.contig_id.isin(filtered_contig), 'clonotype_id'] = ''

    df_all_contig.to_csv(f'{sample}_all_contig.csv', sep=',', index=False)
    df_filter_contig.to_csv(f'{sample}_filtered_contig.csv', sep=',', index=False)
        

def gen_summary(sample, fq2, assembled_reads, seqtype, df_for_clono):
    """ Generate metrics in html 
    """
    df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
    cell_barcodes = set(df_for_clono_pro['barcode'])
    total_cells =len(cell_barcodes)

    read_count = 0
    read_count_all = 0
    umi_dict = defaultdict(set)
    umi_count = defaultdict()
    fq = pyfastx.Fastx(fq2)
    for read in fq:
        read_count_all+=1
        cb = read.name.split(':')[0]
        umi = read.name.split(':')[1]
        umi_dict[cb].add(umi)
        if cb in cell_barcodes:
            read_count+=1
    for cb in umi_dict:
        umi_count[cb] = len(umi_dict[cb])
    df_umi = pd.DataFrame.from_dict(umi_count, orient='index', columns=['UMI'])
    df_umi['barcode'] = df_umi.index
    df_umi = df_umi.reset_index(drop=True)
    df_umi = df_umi.reindex(columns=['barcode', 'UMI'])
    df_umi = df_umi.sort_values(by='UMI', ascending=False)
    df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in cell_barcodes else 'UB')
    df_umi.to_csv(f'{sample}.count.txt', sep='\t', index=False)
    # self.add_data(chart=get_plot_elements.plot_barcode_rank(f'{self.outdir}/count.txt'))

    used_read = 0
    fa = pyfastx.Fastx(assembled_reads)
    for read in fa:
        bc = read.name.split(':')[0]
        if bc in cell_barcodes:
            used_read += 1
    
    if seqtype == 'TCR':
        chains = ['TRA', 'TRB']
    else:
        chains = ['IGH', 'IGL', 'IGK']

    data_dict = {
        "Estimated Number of Cells": total_cells,
        "Mean Read Pairs per Cell": int(read_count/total_cells),
        "Mean Used Read Pairs per Cell": int(used_read/total_cells),
        "Fraction of Reads in Cells": f'({round(used_read / read_count_all * 100, 2)}%)',
    }

    for c in chains:
        temp_df = df_for_clono_pro[df_for_clono_pro['chain']==c]

        try:
            median_umi = int(temp_df['umis'].median())
        except ValueError:
                # ValueError: cannot convert float NaN to integer
            median_umi = 0

        data_dict.update({f'Median {c} UMIs per Cell' : median_umi})


    stats_file = sample + ".cells_stats.json"
    with open(stats_file, "w") as f:
        json.dump(data_dict, f)
        
    # UMI count
    umi_count.loc[lambda x: x > 0]
    umi_count = umi_count.sort_values(ascending=False)
    plot_data = {}
    n = len(umi_count)
    first_noncell = n - 1
    for i, bc in enumerate(umi_count.index):
        if bc not in cell_barcodes:
            first_noncell = i
            break
    last_cell = 0
    for i in range(min(n - 1, MAX_CELL), -1, -1):
        bc = umi_count.index[i]
        if bc in cell_barcodes:
            last_cell = i
            break
    pure = sample + ".cells.pure" + f"({first_noncell}/{first_noncell}, 100%)"
    bg = sample + ".cells.background"
    plot_data[pure] = {}
    plot_data[bg] = {}
    for i in range(first_noncell):
        plot_data[pure][i + 1] = int(umi_count.iloc[i])

    n_mix = last_cell - first_noncell + 1
    if n_mix != 0:
        n_total = len(cell_barcodes)
        n_mix_cell = n_total - first_noncell
        mix_rate = round(n_mix_cell / n_mix * 100, 2)
        mix = sample + ".cells.mix" + f"({n_mix_cell}/{n_mix}, {mix_rate}%)"
        plot_data[mix] = {}
        for i in range(first_noncell, last_cell + 1):
            plot_data[mix][i + 1] = int(umi_count.iloc[i])

    for i in range(last_cell + 1, min(MAX_CELL, n), 10):
        plot_data[bg][i + 1] = int(umi_count.iloc[i])
    # do not record every umi count
    for i in range(MAX_CELL, n, 1000):
        plot_data[bg][i + 1] = int(umi_count.iloc[i])

    umi_file = sample + ".umi_count.json"
    with open(umi_file, "w") as f:
        json.dump(plot_data, f)

   
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--seqtype', required=True)
    parser.add_argument('--coef', required=True) # help='coef for auto filter', default=5
    parser.add_argument("--expected_target_cell_num", required=True)
    #help="Expected T or B cell number. If `--target_cell_barcode` is provided, this argument is ignored.", 
    #type=int,
    #default=3000,
    parser.add_argument('--target_cell_barcode')
    #help="Barcode of target cells. Auto or path of plain text file with one barcode per line",
    #default=None)
    parser.add_argument("--target_weight", required=True)
    #help="UMIs of the target cells are multiplied by this factor. Only used when `--target_cell_barcode` is provided.", 
    #type=float,
    #default=6.0,
    parser.add_argument('--fq2', required=True) #  help='Barcode R2 reads.'
    parser.add_argument('--assembled_reads', required=True)
    parser.add_argument('--filter_report_tsv', required=True)
    parser.add_argument('--annot_fa', required=True)
    parser.add_argument('--barcode_report', required=True)
    args = parser.parse_args()
    
    original_df = parse_contig_file(args.sample, args.barcode_report, args.annot_fa)
    df_for_clono, cell_barcodes, filtered_contig = cell_calling(
        original_df, args.seqtype, args.filter_report_tsv,
        args.expected_target_cell_num, args.target_barcodes, args.target_weight, args.coef
    )
    filter_fasta(args.sample, cell_barcodes)
    parse_clonotypes(args.sample, original_df, df_for_clono, cell_barcodes, filtered_contig)
    gen_summary(args.sample, args.fq2, args.assembled_reads, args.seqtype, df_for_clono)
    
    
