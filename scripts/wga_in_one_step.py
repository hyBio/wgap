#!/home/huyan/miniconda3/bin/python
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
from multiprocessing import Pool
import pandas as pd

_version = "1.0.0"


class wga_in_one_step(object):

    def __init__(self):
        description = 'wga_in_one_step: whole genome alignment in one step with last\n' \
                      'https://github.com/hyBio/L.guttatus/blob/master/scripts/wga_in_one_step.py\n'
        self.parser = argparse.ArgumentParser(
            prog='wga_in_one_step',
            description=description)
        self.parser.add_argument('-b', '--begin', default="fasta_swap", help='the begin subprocess of the pipeline\n'
                                                                             'default: fasta_swap\n'
                                                                             'optional parameters:fasta_download\tfasta_swap\tlastdb\tlast_train\tlastal\tsort\tmultiz\tmaf2lst_fa\tiqtree\n',
                                 required=True)
        self.parser.add_argument('-p', '--parallel', default=1, type=int, help='the number of parallel subprocess')
        self.parser.add_argument('-t', '--threads', type=int, default=1,
                                 help='the number of threads for each subprocess')
        self.parser.add_argument('-c', '--configure', help='configure file path')
        self.parser.add_argument('-a', '--accession_list', help='accession number list file path', default="")
        self.parser.add_argument('-f', '--fa_dir', help='fasta file directory', default=os.getcwd())
        self.parser.add_argument('-r', '--ref_fa', help='reference fasta file')
        self.parser.add_argument('-e', '--exclude_fa', default="", help='exclude fasta file name')
        self.parser.add_argument('-o', '--out_dir', help='output directory', default=os.getcwd())
        self.parser.add_argument('-lastdb', '--lastdb', help='lastdb argparse', type=str, default='')
        self.parser.add_argument('-last_train', '--last_train', help='last_train argparse', type=str, default='')
        self.parser.add_argument('-lastal', '--lastal', help='lastal argparse', type=str, default='')
        self.parser.add_argument('-iqtree', '--iqtree', help='iqtree argparse', type=str, default='')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + _version)

        if len(sys.argv) == 1:
            self.parser.print_help()
            sys.exit(1)
        else:
            self.args = self.parser.parse_args()
            whole_genome_alignment(self.args)
            run(self.args)

class whole_genome_alignment(object):

    def __init__(self, args):
        self.args = args
        self.parallel = self.args.parallel
        self.begin = self.args.begin
        self.threads = self.args.threads
        self.configure = self.args.configure
        self.accession_list = self.args.accession_list
        self.fa_dir = self.args.fa_dir
        self.ref_fa = self.args.ref_fa
        self.exclude_fa = self.args.exclude_fa
        self.out_dir = self.args.out_dir
        self.lastdb = self.args.lastdb
        self.last_train = self.args.last_train
        self.lastal = self.args.lastal
        self.iqtree = self.args.iqtree

        self.acc_list = []
        self.fna_gz_list = []
        self.fasta_swap_name_list = []
        self.exclude_fa_list = []

        if self.accession_list != "":
            with open(self.accession_list) as f:
                lines = f.readlines()
            for line in lines:
                self.acc_list.append(line.strip())

        if self.exclude_fa != "":
            with open(self.exclude_fa, 'r') as f:
                lines = f.readlines()
            for line in lines:
                self.exclude_fa_list.append(line.strip())

        with open(self.configure, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith('#') or line.startswith('\n'):
                continue
            else:
                self.fna_gz_list.append(line.split('\t')[0].strip())
                self.fasta_swap_name_list.append(line.split('\t')[1].strip())
        self.name_list = [i for i in self.fasta_swap_name_list if i != self.ref_fa]
        for i in self.exclude_fa_list:
            self.name_list.remove(i)
        self.last_train_name_list = self.name_list
        self.lastal_name_list = self.name_list
        self.sort_name_list = self.name_list
        self.multiz_name_list = self.name_list

    def fasta_download_func(self, accession):
        self.accession = accession
        with open("{}_logfile:".format(self.accession), 'w', encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directiry:{}\n".format(self.fa_dir))
            f.write(
                "\nrsync --copy-links recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/{}* {}\n".format(
                    self.accession[0:3], self.accession[4:7], self.accession[7:10], self.accession[10:13],
                    self.accession, self.fa_dir, self.accession))
        f.close()

        self.sp_fasta_download = sp.Popen(shlex.split(
            "rsync --copy-links recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/{}* {}".format(
                self.accession[0:3], self.accession[4:7], self.accession[7:10], self.accession[10:13], self.accession,
                self.fa_dir)), stdout=sp.PIPE, stderr=sp.PIPE)
        self.sp_fasta_download.communicate()

        if self.sp_fasta_download.returncode == 0:
            with open("{}_logfile:".format(self.accession), 'a', encoding='utf-8') as f:
                f.write(self.sp_fasta_download.stdout.read().decode('utf-8'))
                f.write("\nfasta_download finished\n")
            f.close()
            return 0
        else:
            with open("error_{}:".format(self.accession), 'a', encoding='utf-8') as f:
                f.write("\nfasta_download failed\n" + "\nfasta download failed message:\n")
                f.write(self.sp_fasta_download.stderr.read().decode('utf-8'))
            f.close()
            return 1

    def fasta_swap_func(self, fna_gz_fasta_swap_name):
        fna_gz = fna_gz_fasta_swap_name[0]
        fasta_swap_name = fna_gz_fasta_swap_name[1]

        sp.Popen(shlex.split("mkdir -p {}/00_assembly_fasta".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/00_assembly_fasta/".format(self.out_dir)))
        os.chdir(r"{}00_assembly_fasta/".format(self.out_dir))

        with open("{}/00_assembly_fasta/{}_logfile".format(self.out_dir, fasta_swap_name), "w",
                  encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write(" ".join(
                ['\nrunning command line:', 'zcat', '{}'.format(fna_gz), '|', 'awk', '\'{print $1}\'', '|', 'sed',
                 '\'s/>/>{}_/g\''.format(fasta_swap_name), '>',
                 '{}/00_assembly_fasta/{}.fa\n'.format(self.out_dir, fasta_swap_name)]))
        f.close()

        self.sp_ungzip = sp.Popen(["zcat", "{}".format(fna_gz)], stdout=sp.PIPE)
        self.sp_awk = sp.Popen(shlex.split("awk '{print $1}'"), stdin=self.sp_ungzip.stdout, stdout=sp.PIPE)
        self.sp_ungzip.stdout.close()
        self.sp_sed = sp.Popen(shlex.split("sed 's/>/>{}_/g'".format(fasta_swap_name)), stdin=self.sp_awk.stdout,
                               stdout=sp.PIPE)
        self.sp_awk.stdout.close()
        fasta_swap_1_output = self.sp_sed.communicate()[0]
        with open('{}/00_assembly_fasta/{}.fa'.format(self.out_dir, fasta_swap_name), 'w', encoding='utf-8') as fa:
            fa.write(fasta_swap_1_output.decode('UTF-8').strip() + '\n')
        self.sp_sed.stdout.close()

        if self.sp_sed.returncode == 0:
            with open("{}/00_assembly_fasta/{}_logfile".format(self.out_dir, fasta_swap_name), "a",
                      encoding='utf-8') as f:
                f.write("\nfasta_swap successfully finished\n")
                f.write("\nrunning command line:samtools faidx {}/00_assembly_fasta/{}.fa\n".format(self.out_dir,
                                                                                                    fasta_swap_name))
            f.close()
            self.sp_faidx = sp.Popen(
                shlex.split("samtools faidx {}/00_assembly_fasta/{}.fa".format(self.out_dir, fasta_swap_name)),
                stdout=sp.PIPE)
            fasta_swap_2_output = self.sp_faidx.communicate()[0]
            with open('{}/00_assembly_fasta/{}.fa.fai'.format(self.out_dir, fasta_swap_name), 'w',
                      encoding='utf-8') as f:
                f.write(fasta_swap_2_output.decode('UTF-8').strip() + '\n')
            self.sp_faidx.stdout.close()

            if self.sp_faidx.returncode == 0:
                with open("{}/00_assembly_fasta/{}_logfile".format(self.out_dir, fasta_swap_name), "a",
                          encoding='utf-8') as f:
                    f.write("\nsamtools faidx successfully finished\n")
                f.close()
                return 0
            else:
                with open("{}/00_assembly_fasta/error_{}".format(self.out_dir, fasta_swap_name), "w",
                          encoding='utf-8') as f:
                    f.write("\nsamtools faidx failed\n" + "samtools faidx failed message:\n")
                    f.write(self.sp_faidx.stderr.read().decode('UTF-8'))
                f.close()
                return 1
        else:
            with open("{}/00_assembly_fasta/error_{}".format(self.out_dir, fasta_swap_name), "w",
                      encoding='utf-8') as f:
                f.write("\nsed failed\n" + "sed failed message:\n")
                f.write(self.sp_sed.stderr.read().decode('UTF-8'))
            f.close()
            return 1

    def lastdb_func(self):
        sp.Popen(shlex.split("mkdir -p {}/01_lastdb".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/01_lastdb/".format(self.out_dir)))
        os.chdir(r"{}/01_lastdb/".format(self.out_dir))

        with open("{}/01_lastdb/logfile".format(self.out_dir), "w", encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write("\nrunning command line:lastdb -P {} ".format(
                self.parallel * self.threads) + self.lastdb + " " + self.ref_fa + '_db ' + self.out_dir + '/00_assembly_fasta/' + self.ref_fa + '.fa\n')
        f.close()

        sp_lastdb_cmd = "lastdb -P {} ".format(
            self.parallel * self.threads) + self.lastdb + " " + self.ref_fa + '_db ' + self.out_dir + '/00_assembly_fasta/' + self.ref_fa + '.fa'
        self.sp_lastdb = sp.Popen(shlex.split(sp_lastdb_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
        # wiat for the process to finish
        self.sp_lastdb.communicate()
        # check lastdb return code
        if self.sp_lastdb.returncode == 0:
            with open("{}/01_lastdb/logfile".format(self.out_dir), "a", encoding='utf-8') as f:
                f.write("\nlastdb successfully finished\n")
            f.close()
            return 0
        else:
            self.sp_lastdb_help = sp.Popen(shlex.split("lastdb --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            with open("{}/01_lastdb/logfile".format(self.out_dir), "a", encoding='utf-8') as f:
                f.write("\nlastdb failed\n" + "lastdb failed message:\n")
                f.write(self.sp_lastdb.stderr.read().decode('UTF-8'))
                f.write("\nlastdb help:\n")
                f.write(self.sp_lastdb_help.stdout.read().decode('UTF-8'))
                f.write(self.sp_lastdb_help.stderr.read().decode('UTF-8'))
            f.close()
            return 1

    def last_train_func(self, last_train_name):
        sp.Popen(shlex.split("mkdir -p {}/02_last_train".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/02_last_train/".format(self.out_dir)))
        os.chdir(r"{}/02_last_train/".format(self.out_dir))

        with open("{}/02_last_train/{}_logfile".format(self.out_dir, last_train_name), "w", encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write(
                "\nrunning command line:last-train -P {} {} {}/01_lastdb/{}_db {}/00_assembly_fasta/{}.fa > {}/02_last_train/{}.mat\n".format(
                    self.threads, self.last_train, self.out_dir, self.ref_fa, self.out_dir, last_train_name,
                    self.out_dir, last_train_name))
        f.close()

        sp_last_train_cmd = "last-train -P {} ".format(
            self.threads) + self.last_train + " " + self.out_dir + "/01_lastdb/" + self.ref_fa + "_db " + self.out_dir + "/00_assembly_fasta/" + last_train_name + ".fa"
        self.sp_last_train = sp.Popen(shlex.split(sp_last_train_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
        # wiat for the process to finish
        last_train_output = self.sp_last_train.communicate()[0]
        with open("{}/02_last_train/{}.mat".format(self.out_dir, last_train_name), "w") as f:
            f.write(last_train_output.decode('UTF-8').strip() + "\n")
        self.sp_last_train.stdout.close()

        # check lastdb return code
        if self.sp_last_train.returncode == 0:
            with open("{}/02_last_train/{}_logfile".format(self.out_dir, last_train_name), "a",
                      encoding='utf-8') as f:
                f.write("\nlast-train successfully finished\n")
            f.close()
            return 0
        else:
            with open("{}/02_last_train/error_{}".format(self.out_dir, last_train_name), "w",encoding='utf-8') as f:
                f.write("\nlast-train failed\n" + "last-train error message:\n")
                f.write(self.sp_last_train.stderr.read().decode('UTF-8'))
            f.close()
            self.sp_last_train.stderr.close()
            return 1

    def lastal_func(self, lastal_name):
        sp.Popen(shlex.split("mkdir -p {}/03_lastal".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/03_lastal/".format(self.out_dir)))
        os.chdir(r"{}/03_lastal/".format(self.out_dir))

        with open("{}/03_lastal/{}_logfile".format(self.out_dir, lastal_name), "w") as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write(
                "\nrunning command line:lastal -P {} {} {}/02_last_train/{}.mat {}/01_lastdb/{}_db {}/00_assembly_fasta/{}.fa | last-split -f MAF+ > {}/03_lastal/{}.maf\n".format(
                    self.threads, self.lastal, self.out_dir, lastal_name, self.out_dir, self.ref_fa, self.out_dir,
                    lastal_name, self.out_dir, lastal_name))
        f.close()

        sp_lastal_1_cmd = "lastal -P {} ".format(
            self.threads) + self.lastal + " " + self.out_dir + "/02_last_train/" + lastal_name + ".mat " + self.out_dir + "/01_lastdb/" + self.ref_fa + "_db " + self.out_dir + "/00_assembly_fasta/" + lastal_name + ".fa "
        sp_lastal_2_cmd = "last-split " + "-f " + "MAF+"

        self.sp_lastal_1 = sp.Popen(shlex.split(sp_lastal_1_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
        self.sp_lastal_2 = sp.Popen(shlex.split(sp_lastal_2_cmd), stdin=self.sp_lastal_1.stdout, stdout=sp.PIPE,
                                    stderr=sp.PIPE)
        sp_lastal_2_output = self.sp_lastal_2.communicate()[0]
        with open("{}/03_lastal/{}.maf".format(self.out_dir, lastal_name), "w", encoding='utf-8') as f:
            f.write(sp_lastal_2_output.decode('UTF-8').strip() + "\n")
        f.close()
        self.sp_lastal_1.stdout.close()

        # check lastal_2 return code
        if self.sp_lastal_2.returncode == 0:
            with open("{}/03_lastal/{}_logfile".format(self.out_dir, lastal_name), "a", encoding='utf-8') as f:
                f.write("\nlastal and last-split successfully finished\n")
            f.close()
            return 0
        else:
            with open("{}/03_lastal/error_{}".format(self.out_dir, lastal_name), "w", encoding='utf-8') as f:
                f.write("\nerror:lastal and last-split failed\n" + "lastal error message:\n")
                f.write(self.sp_lastal_1.stderr.read().decode('UTF-8'))
                f.write("\nlast-split error message:\n")
                f.write(self.sp_lastal_2.stderr.read().decode('UTF-8'))
            f.close()
            self.sp_lastal_1.stderr.close()
            self.sp_lastal_2.stderr.close()
            return 1

    def sort_func(self, sort_name):
        sp.Popen(shlex.split("mkdir -p {}/04_sort".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/04_sort/".format(self.out_dir)))
        os.chdir(r"{}/04_sort/".format(self.out_dir))

        with open("{}/04_sort/{}_logfile".format(self.out_dir, sort_name), "w") as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write(
                "\nrunning command line:maf-swap {}/03_lastal/{}.maf |last-split |maf-swap |maf-sort|grep -v '#'|sed '1i\#maf version=1.0 scoring=last' > {}/04_sort/{}.maf\n".format(
                    self.out_dir, sort_name, self.out_dir, sort_name))
        f.close()

        sp_sort_1_cmd = "maf-swap {}/03_lastal/{}.maf".format(self.out_dir, sort_name)
        sp_sort_2_cmd = "last-split"
        sp_sort_3_cmd = "maf-swap"
        sp_sort_4_cmd = "maf-sort"

        self.sp_sort_1 = sp.Popen(shlex.split(sp_sort_1_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
        self.sp_sort_2 = sp.Popen(shlex.split(sp_sort_2_cmd), stdin=self.sp_sort_1.stdout, stdout=sp.PIPE,
                                  stderr=sp.PIPE)
        self.sp_sort_3 = sp.Popen(shlex.split(sp_sort_3_cmd), stdin=self.sp_sort_2.stdout, stdout=sp.PIPE,
                                  stderr=sp.PIPE)
        self.sp_sort_4 = sp.Popen(shlex.split(sp_sort_4_cmd), stdin=self.sp_sort_3.stdout, stdout=sp.PIPE,
                                  stderr=sp.PIPE)
        sp_sort_4_output = self.sp_sort_4.communicate()[0]
        with open("{}/04_sort/{}.maf".format(self.out_dir, sort_name), "w", encoding='utf-8') as f:
            f.write("##maf version=1 scoring=last\n")
            for line in sp_sort_4_output.decode('UTF-8').split('\n'):
                if not line.startswith("#"):
                    f.write(line + "\n")
        f.close()
        self.sp_sort_1.stdout.close()
        self.sp_sort_2.stdout.close()
        self.sp_sort_3.stdout.close()
        self.sp_sort_4.stdout.close()

        # check sort_4 return code
        if self.sp_sort_4.returncode == 0:
            with open("{}/04_sort/{}_logfile".format(self.out_dir, sort_name), "a", encoding='utf-8') as f:
                f.write("\nmaf-sort successfully finished\n")
            f.close()
            return 0
        else:
            with open("{}/04_sort/error_{}".format(self.out_dir, sort_name), "w", encoding='utf-8') as f:
                f.write("\nerror:maf-sort failed\n" + "maf-swap error message:\n")
                f.write(self.sp_sort_1.stderr.read().decode('UTF-8'))
                f.write("\nlast-split error message:\n")
                f.write(self.sp_sort_2.stderr.read().decode('UTF-8'))
                f.write("\nmaf-swap error message:\n")
                f.write(self.sp_sort_3.stderr.read().decode('UTF-8'))
                f.write("\nmaf-sort error message:\n")
                f.write(self.sp_sort_4.stderr.read().decode('UTF-8'))
            f.close()
            self.sp_sort_4.stderr.close()
            self.sp_sort_3.stderr.close()
            self.sp_sort_2.stderr.close()
            self.sp_sort_1.stderr.close()
            return 1

    def multiz_func(self, multiz_name):
        first_name = multiz_name[0]
        sp.Popen(shlex.split("mkdir -p {}/05_multiz".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/05_multiz/".format(self.out_dir)))
        os.chdir(r"{}/05_multiz/".format(self.out_dir))
        sp.Popen(shlex.split(
            "cp {}/04_sort/{}.maf {}/05_multiz/1.maf".format(self.out_dir, first_name, self.out_dir))).wait()

        self.name_n = 1
        for next_name in multiz_name[1::]:
            with open("{}/05_multiz/logfile".format(self.out_dir), "w") as f:
                f.write("\nargparse:{}\n".format(self.args))
                f.write("\nrunning directory:{}\n".format(self.out_dir))
                f.write(
                    "\nrunning command line:multiz {}/05_multiz/{}.maf {}/04_sort/{}.maf 0 U1 U2 > {}/05_multiz/{}.maf\n".format(self.out_dir,self.name_n,self.out_dir,next_name,self.out_dir,self.name_n+1))
            f.close()

            sp_multiz_cmd = "multiz {}/05_multiz/{}.maf {}/04_sort/{}.maf 0 U1 U2".format(self.out_dir,self.name_n,self.out_dir,next_name)
            self.sp_multiz = sp.Popen(shlex.split(sp_multiz_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
            sp_multiz_output = self.sp_multiz.communicate()[0]
            with open("{}/05_multiz/{}.maf".format(self.out_dir, self.name_n+1), "w", encoding='utf-8') as f:
                f.write(sp_multiz_output.decode('UTF-8'))
            f.close()

            if self.sp_multiz.returncode == 0:
                with open("{}/05_multiz/logfile".format(self.out_dir), "a", encoding='utf-8') as f:
                    f.write("\n{} multiz successfully finished\n".format(self.name_n+1))
                f.close()
            else:
                with open("{}/05_multiz/error_log".format(self.out_dir), "w", encoding='utf-8') as f:
                    f.write("\nerror:multiz failed\n" + "multiz error message:\n")
                    f.write(self.sp_multiz.stderr.read().decode('UTF-8'))
                f.close()
                self.sp_multiz.stderr.close()
                return 1
            self.name_n += 1
        return 0

    def maf2lst_fa_func(self,maf_file,lst_file,fasta_file):
        sp.Popen(shlex.split("mkdir -p {}/06_maf2lst_fa".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/06_maf2lst_fa/".format(self.out_dir)))
        os.chdir(r"{}/06_maf2lst_fa/".format(self.out_dir))

        csplit_cmd = 'csplit %s /^a/ -n8 -s {*} -f block_' %(maf_file)
        sp_csplit = sp.Popen(shlex.split(csplit_cmd))
        sp_csplit.wait()

        for maf_block_file in [f for f in os.listdir('.') if f.startswith('block')]:
            with open(maf_block_file, 'r', encoding='utf-8') as maf_block:
                maf_block_startswith_s = [line for line in maf_block.readlines() if line.startswith('s')]

            if len(maf_block_startswith_s) == self.name_n:
                with open('full.maf', 'w', encoding='utf-8') as full_maf, open(maf_block_file, 'r', encoding='utf-8') as maf_block:
                    full_maf.write(maf_block.read())
            else:
                with open('broken.maf', 'w', encoding='utf-8') as broken_maf, open(maf_block_file, 'r', encoding='utf-8') as maf_block:
                    broken_maf.write(maf_block.read())
            os.remove(maf_block_file)

        with open("full.maf", 'r') as maf_f, open(fasta_file, 'w') as fa_f:
            maf = [' '.join(line.split()).split() for line in maf_f.readlines() if line.startswith('s')]
            maf_df = pd.DataFrame(maf)
            maf_df.columns = ['s', 'chr_name', 'start', 'base_len', 'strand', 'chr_len', 'seq']
            maf_df['species'] = maf_df.chr_name.apply(lambda x: x.split('_')[0])

            lst_df = maf_df.loc[maf_df.species == maf_df.species[0], :]
            lst_df = lst_df.apply(lambda x:pd.DataFrame({"chr":[x.chr_name for i in range(len(x.seq))],"start":[int(x.start)+i for i in range(len(x.seq))]}),axis=1)
            lst_df = pd.concat([i for i in lst_df.values], axis=0)
            for species in maf_df.species.unique():
                species_seq = [i for i in ''.join(maf_df.loc[maf_df.species == species, 'seq'].values).upper()]
                lst_df[species] = species_seq
            lst_df.to_csv(lst_file, sep='\t', header=True, index=False)

            species_seq = maf_df.groupby("species")['seq'].apply(lambda x: ''.join(x))
            for species, seq in species_seq.items():
                fa_f.write('>' + species + '\n')
                fa_f.write(seq + '\n')
        return 0

    def iqtree_func(self,fasta_file):
        sp.Popen(shlex.split("cd {}/06_maf2lst_fa/".format(self.out_dir)))
        os.chdir(r"{}/06_maf2lst_fa/".format(self.out_dir))

        with open("iqtree_logfile", "w", encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write("\nrunning command line:iqtree -s {} -nt {} {}\n".format(fasta_file, self.threads, self.iqtree))
        f.close()

        self.sp_iqtree = sp.Popen(shlex.split("iqtree -s {} -nt {} {}".format(fasta_file, self.threads, self.iqtree)), stdout=sp.PIPE, stderr=sp.PIPE)
        self.sp_iqtree.wait()

        if self.sp_iqtree.returncode == 0:
            with open("iqtree_logfile", "a", encoding='utf-8') as f:
                f.write("\n{} iqtree successfully finished\n".format(self.name_n))
                f.write(self.sp_iqtree.stdout.read().decode('UTF-8'))
            f.close()
            return 0
        else:
            with open("iqtree_error_log", "w", encoding='utf-8') as f:
                f.write("\nerror:iqtree failed\n" + "iqtree error message:\n")
                f.write(self.sp_iqtree.stderr.read().decode('UTF-8'))
            f.close()
            self.sp_iqtree.stderr.close()
            return 1


class run(whole_genome_alignment):

    def __init__(self, args):
        super(run, self).__init__(args)
        if self.begin == "fasta_download":
            self.fasta_download_run()
        elif self.begin == "fasta_swap":
            if self.fasta_swap_run() == 0:
                if self.lastdb_run() == 0:
                    if self.last_train_run() == 0:
                        if self.lastal_run() == 0:
                            if self.sort_run() == 0:
                                if self.multiz_run() == 0:
                                    if self.maf2lst_fa_run() == 0:
                                        if self.iqtree_run() == 0:
                                            print("All done!")
                                            sys.exit(0)
                                        else:
                                            print("iqtree finished with error!")
                                            sys.exit(1)
                                    else:
                                        print("maf2lst_fa finished with error!")
                                        sys.exit(1)
                                else:
                                    print("multiz_run finished with error!")
                                    sys.exit(1)
                            else:
                                print("sort finished with error!")
                                sys.exit(1)
                        else:
                            print("lastal finished with error!")
                            sys.exit(1)
                    else:
                        print("last_train finished with error!")
                        sys.exit(1)
                else:
                    print("lastdb finished with error!")
                    sys.exit(1)
            else:
                print("fasta_swap finished with error!")
                sys.exit(1)
        elif self.begin == "lastdb":
            if self.lastdb_run() == 0:
                if self.last_train_run() == 0:
                    if self.lastal_run() == 0:
                        if self.sort_run() == 0:
                            if self.multiz_run() == 0:
                                if self.maf2lst_fa_run() == 0:
                                    if self.iqtree_run() == 0:
                                        print("All done!")
                                        sys.exit(0)
                                    else:
                                        print("iqtree finished with error!")
                                        sys.exit(1)
                                else:
                                    print("maf2lst_fa finished with error!")
                                    sys.exit(1)
                            else:
                                print("multiz_run finished with error!")
                                sys.exit(1)
                        else:
                            print("sort finished with error!")
                            sys.exit(1)
                    else:
                        print("lastal finished with error!")
                        sys.exit(1)
                else:
                    print("last_train finished with error!")
                    sys.exit(1)
            else:
                print("lastdb finished with error!")
                sys.exit(1)
        elif self.begin == "last_train":
            if self.last_train_run() == 0:
                if self.lastal_run() == 0:
                    if self.sort_run() == 0:
                        if self.multiz_run() == 0:
                            if self.maf2lst_fa_run() == 0:
                                if self.iqtree_run() == 0:
                                    print("All done!")
                                    sys.exit(0)
                                else:
                                    print("iqtree finished with error!")
                                    sys.exit(1)
                            else:
                                print("maf2lst_fa finished with error!")
                                sys.exit(1)
                        else:
                            print("multiz_run finished with error!")
                            sys.exit(1)
                    else:
                        print("sort finished with error!")
                        sys.exit(1)
                else:
                    print("lastal finished with error!")
                    sys.exit(1)
            else:
                print("last_train finished with error!")
                sys.exit(1)
        elif self.begin == "lastal":
            if self.lastal_run() == 0:
                if self.sort_run() == 0:
                    if self.multiz_run() == 0:
                        if self.maf2lst_fa_run() == 0:
                            if self.iqtree_run() == 0:
                                print("All done!")
                                sys.exit(0)
                            else:
                                print("iqtree finished with error!")
                                sys.exit(1)
                        else:
                            print("maf2lst_fa finished with error!")
                            sys.exit(1)
                    else:
                        print("multiz_run finished with error!")
                        sys.exit(1)
                else:
                    print("sort finished with error!")
                    sys.exit(1)
            else:
                print("lastal finished with error!")
                sys.exit(1)
        elif self.begin == "sort":
            if self.sort_run() == 0:
                if self.multiz_run() == 0:
                    if self.maf2lst_fa_run() == 0:
                        if self.iqtree_run() == 0:
                            print("All done!")
                            sys.exit(0)
                        else:
                            print("iqtree finished with error!")
                            sys.exit(1)
                    else:
                        print("maf2lst_fa finished with error!")
                        sys.exit(1)
                else:
                    print("multiz_run finished with error!")
                    sys.exit(1)
            else:
                print("sort finished with error!")
                sys.exit(1)
        elif self.begin == "multiz":
            if self.multiz_run() == 0:
                if self.maf2lst_fa_run() == 0:
                    if self.iqtree_run() == 0:
                        print("All done!")
                        sys.exit(0)
                    else:
                        print("iqtree finished with error!")
                        sys.exit(1)
                else:
                    print("maf2lst_fa finished with error!")
                    sys.exit(1)
            else:
                print("multiz_run finished with error!")
                sys.exit(1)
        elif self.begin == "maf2lst_fa":
            if self.maf2lst_fa_run() == 0:
                if self.iqtree_run() == 0:
                    print("All done!")
                    sys.exit(0)
                else:
                    print("iqtree finished with error!")
                    sys.exit(1)
            else:
                print("maf2lst_fa finished with error!")
                sys.exit(1)
        else:
            print("Incorrect parameter '{}' given to whole_genome_alignmen, please check the help".format(self.begin))
            sys.exit(1)

    def fasta_download_run(self):
        p = Pool(self.parallel)
        result = []
        for i in self.acc_list:
            result.append(p.apply_async(self.fasta_download_func, args=(i,)))
        p.close()
        p.join()
        [i.wait() for i in result]
        if len([i.get() for i in result if i.ready() and i.successful()]) == len(self.acc_list):
            if sum([i.get() for i in result if i.ready() and i.successful()]) == 0:
                return 0
            else:
                return 1

    def fasta_swap_run(self):
        p = Pool(self.parallel)
        result = []
        for i in [(x,y) for x,y in zip(self.fna_gz_list, self.fasta_swap_name_list)]:
            result.append(p.apply_async(self.fasta_swap_func, args=(i,)))
        p.close()
        p.join()
        [i.wait() for i in result]
        if len([i.get() for i in result if i.ready() and i.successful()]) == len(self.fna_gz_list):
            if sum([i.get() for i in result if i.ready() and i.successful()]) == 0:
                return 0
            else:
                return 1

    def lastdb_run(self):
        if self.lastdb_func() == 0:
            return 0
        else:
            return 1

    def last_train_run(self):
        p = Pool(self.parallel)
        result = []
        for i in self.last_train_name_list:
            result.append(p.apply_async(self.last_train_func, args=(i,)))
        p.close()
        p.join()
        [i.wait() for i in result]
        print([i.get() for i in result])
        if len([i.get() for i in result if i.ready() and i.successful()]) == len(self.last_train_name_list):
            if sum([i.get() for i in result if i.ready() and i.successful()]) == 0:
                return 0
            else:
                return 1

    def lastal_run(self):
        p = Pool(self.parallel)
        result = []
        for i in self.lastal_name_list:
            result.append(p.apply_async(self.lastal_func, args=(i,)))
        p.close()
        p.join()
        [i.wait() for i in result]
        if len([i.get() for i in result if i.ready() and i.successful()]) == len(self.lastal_name_list):
            if sum([i.get() for i in result if i.ready() and i.successful()]) == 0:
                return 0
            else:
                return 1

    def sort_run(self):
        p = Pool(self.parallel)
        result = []
        for i in self.sort_name_list:
            result.append(p.apply_async(self.sort_func, args=(i,)))
        p.close()
        p.join()
        [i.wait() for i in result]
        if len([i.get() for i in result if i.ready() and i.successful()]) == len(self.sort_name_list):
            if sum([i.get() for i in result if i.ready() and i.successful()]) == 0:
                return 0
            else:
                return 1

    def multiz_run(self):
        if self.multiz_func(self.multiz_name_list) == 0:
            return 0
        else:
            return 1

    def maf2lst_fa_run(self):
        if not hasattr(self,'name_n'):
            name_n = len(self.fna_gz_list)
        else:
            name_n = self.name_n
        maf_file = "{}/05_multiz/{}.maf".format(self.out_dir, name_n-1)
        lst_file = "{}/06_maf2lst_fa/{}_{}.lst".format(self.out_dir, self.ref_fa, name_n-1)
        fasta_file = "{}/06_maf2lst_fa/{}_{}.fa".format(self.out_dir, self.ref_fa, name_n-1)
        if self.maf2lst_fa_func(maf_file, lst_file, fasta_file) == 0:
            return 0
        else:
            return 1

    def iqtree_run(self):
        if not hasattr(self,'name_n'):
            name_n = len(self.fna_gz_list)
        else:
            name_n = self.name_n
        fasta_file = "{}/06_maf2lst_fa/{}_{}.fa".format(self.out_dir, self.ref_fa, name_n-1)
        if self.iqtree_func(fasta_file) == 0:
            return 0
        else:
            return 1


if __name__ == '__main__':
    wga_in_one_step()



    