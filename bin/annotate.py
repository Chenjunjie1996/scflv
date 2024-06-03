#!/usr/bin/env python

import pandas as pd
import argparse
import json
import utils

    
def get_vdj_metric(df, chains, pairs):
    """
    Add vdj metrics.
    """
    data_dict = {}
    fl_pro_pair_df = pd.DataFrame(df[df['productive']==True].barcode.value_counts())
    fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode']>=2]
    cell_nums = len(set(df['barcode']))
    
    data_dict = {
        'Cells With Productive V-J Spanning Pair': utils.format_value(fl_pro_pair_df.shape[0], cell_nums)
    }

    for pair in pairs:
        chain1, chain2 = pair.split('_')[0], pair.split('_')[1]
        cbs1 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain1)].barcode)
        cbs2 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain2)].barcode)
        paired_cbs = len(cbs1.intersection(cbs2))

        data_dict.update({
            'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair': utils.format_value(paired_cbs, cell_nums)
        })

    for chain in chains:
        value = len(set(df[df['chain']==chain].barcode))
        data_dict.update({
            f'Cells With {chain} Contig': utils.format_value(value, cell_nums)
        })

        value = len(set(df[(df['chain']==chain)&(df['cdr3']!=None)].barcode))
        data_dict.update({
            f'Cells With CDR3-annotated {chain} Contig': utils.format_value(value, cell_nums)
        })
        
        value = len(set(df[(df['full_length']==True)&(df['chain']==chain)].barcode))
        data_dict.update({
            f'Cells With V-J Spanning {chain} Contig': utils.format_value(value, cell_nums)
        })

        value = len(set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain)].barcode))
        data_dict.update({
            f'Cells With Productive {chain} Contig': utils.format_value(value, cell_nums)
        })

    return data_dict


def parse_clonotype(df):
    """Generate clonotypes table.
    """
    df_clonotypes=pd.read_csv(df, sep=',')
    df_clonotypes['ClonotypeID'] = df_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
    df_clonotypes['Frequency'] = df_clonotypes['frequency']
    df_clonotypes['Proportion'] = df_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
    df_clonotypes['CDR3_aa'] = df_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))

    return df_clonotypes

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--seqtype', required=True)
    parser.add_argument('--contig_file', required=True)
    parser.add_argument('--clonotype_file', required=True)
    args = parser.parse_args()

    if args.seqtype == 'BCR':
        chains = ['IGH', 'IGL', 'IGK']
        paired_groups = ['IGK_IGH', 'IGL_IGH']
    else:
        chains = ['TRA', 'TRB']
        paired_groups = ['TRA_TRB']
        
    df = pd.read_csv(args.contig_file)
    data_dict = get_vdj_metric(args.sample, df, chains, paired_groups)
    stats_file = args.sample + ".V(D)J_Annotation.json"
    with open(stats_file, "w") as f:
        json.dump(data_dict, f)
        
 
    # plot_data = {}
    # plot_data['df_clonotype'] = {}
    # plot_data['barplot'] = {}
        
    # df_clonotypes = parse_clonotype(args.clonotype_file)
    # df_table = df_clonotypes[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']]
    # table_dict = {}
    # table_dict['title'] = 'Clonetypes'
    # table_dict['table'] = df_table.to_html(
    #     escape=False,
    #     index=False,
    #     table_id='clonetypes',
    #     justify="center")
    # table_dict['id'] = 'clonetypes'
    # plot_data['df_clonotype'] = table_dict

    # df_clonotypes['ClonotypeID'] = df_clonotypes['ClonotypeID'].astype("int")
    # df_clonotypes.sort_values(by=['ClonotypeID'], inplace=True)
    # Barplot = Bar_plot(df_bar=df_clonotypes).get_plotly_div()
    # self.add_data(Barplot=Barplot)