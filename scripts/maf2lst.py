#!/home/huyan/miniconda3/bin/python
# -*- coding: utf-8 -*-
# @Time : 2022/4/16 16:22
# @Author : huyan
# @FileName: maf2lst.py
# @Software: wgap
# @Site :


import sys
import pandas as pd

maf_file = sys.argv[1]
lst_file = sys.argv[2]
fasta_file = sys.argv[3]

class Maf2Lst(object):
    def __init__(self, maf_file, lst_file, fasta_file):
        self.maf_file = maf_file
        self.lst_file = lst_file
        self.fasta_file = fasta_file

    def maf2lst_fa_func(self):
        with open(maf_file, 'r') as maf_f, open(lst_file, 'w') as lst_f, open(fasta_file, 'w') as fa_f:
            maf = [' '.join(line.split()).split() for line in maf_f.readlines() if line.startswith('s')]
            maf_df = pd.DataFrame(maf)
            maf_df.columns = ['s', 'chr_name', 'start', 'base_len', 'strand', 'chr_len', 'seq']
            maf_df['species'] = maf_df.chr_name.apply(lambda x: x.split('_')[0])

            lst_df = maf_df.loc[maf_df.species == maf_df.species[0], :]
            lst_df = lst_df.apply(lambda x:pd.DataFrame({"chr":[x.chr_name for i in range(len(x.seq))],"start":[int(x.start)+i for i in range(len(x.seq))]}),axis=1)
            lst_df = pd.concat([i for i in lst_df.values], axis=0)
            for species in maf_df.species.unique()[1:10]:
                lst_df[species] = [i for i in ''.join(maf_df.loc[maf_df.species == species, 'seq'].values).strip()]

            lst_df.to_csv(lst_f, sep='\t', header=True, index=False)

            species_seq = maf_df.groupby("species")['seq'].apply(lambda x: ''.join(x))
            for species, seq in species_seq.items():
                fa_f.write('>' + species + '\n')
                fa_f.write(seq + '\n')
            fa_f.close()
            lst_f.close()
            maf_f.close()


if __name__ == '__main__':
    Maf2Lst(maf_file, lst_file, fasta_file).maf2lst_fa_func()


