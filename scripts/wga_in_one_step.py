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

_version = "1.0.0"

class wga_in_one_step(Object):
    def __init__(self):
        description = 'wga_in_one_step: whole genome alignment in one step with last'
        self.parser = argparse.ArgumentParser(
            prog='wga_in_one_step({})'.format(_version),
            description=description)
        self.parser.add_argument('-f', '--fa_dir', help='fasta file directory')
        self.parser.add_argument('-c', '--configure', help='configure file path')
        self.parser.add_argument('-o', '--out_dir', help='output directory')
        self.parser.add_argument('-s1', '--lastdb', help='lastdb argparse')
        self.parser.add_argument('-s2', '--last_train', help='last_train argparse')
        self.parser.add_argument('-s3', '--lastal', help='lastal argparse')
        self.parser.add_argument('-s4', '--multiz', help='multiz argparse')



