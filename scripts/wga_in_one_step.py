# _*_ coding: utf-8 _*_
# @Time : 2022/4/10 21:12 
# @Author : 胡琰 
# @Version：V 0.1
# @File : wga_in_one_step.py
# @Site :


import argparse
import sys
import os
import shlex
import subprocess as sp

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import last_wga as lw

_version = "1.0.0"

class wga_in_one_step(object):
    def __init__(self):
        description = 'wga_in_one_step{}: whole genome alignment in one step with last'.format(_version)
        self.parser = argparse.ArgumentParser(
            prog='wga_in_one_step',
            description=description)
        self.parser.add_argument('-t', '--threads', type=int, default=1, help='number of threads')
        self.parser.add_argument('-c', '--configure', help='configure file path', required=True)
        self.parser.add_argument('-f', '--fa_dir', help='fasta file directory',default=os.getcwd())
        self.parser.add_argument('-e', '--exclude_fa', help='exclude fasta file')
        self.parser.add_argument('-o', '--out_dir', help='output directory', default=os.getcwd())
        self.parser.add_argument('-s1', '--lastdb', help='lastdb argparse',type=str,default='')
        self.parser.add_argument('-s2', '--last_train', help='last_train argparse',type=str,default='')
        self.parser.add_argument('-s3', '--lastal', help='lastal argparse',type=str,default='')
        #self.parser.add_argument('-s4', '--multiz', help='multiz argparse',type=str,default='')

        if len(sys.argv) == 1:
            self.parser.print_help()
            sys.exit(1)
        else:
            self.args = self.parser.parse_args()
            self.whole_genome_alignment = whole_genome_alignment(self.args)




class whole_genome_alignment(object):

    def __init__(self,args):
        self.args = args
        self.threads = self.args.threads
        self.configure = self.args.configure
        self.fa_dir = self.args.fa_dir
        self.exclude_fa = self.args.exclude_fa
        self.out_dir = self.args.out_dir
        self.lastdb = self.args.lastdb
        self.last_train = self.args.last_train
        self.lastal = self.args.lastal





if __name__ == '__main__':
    wga_in_one_step()