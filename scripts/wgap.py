#!python
# -*- coding: utf-8 -*-
# @Time : 2022/4/19 23:11
# @Author : huyan
# @FileName: wgap.py
# @Software: wgap
# @Site :


import argparse
import sys
import os
import shlex
import subprocess as sp
import multiprocessing as mp
import pandas as pd

_version = "1.0.0"



class Wgap_help(object):

    def __init__(self):
        description = 'Wgap: whole genome alignment pipe with last\nhttps://github.com/hyBio/L.guttatus/blob/master/scripts/wga_in_one_step.py\n'
        self.parser = argparse.ArgumentParser(prog='Wgap', description=description)
        self.parser.add_argument('-b', '--begin', default="fasta_swap", help='the begin subprocess of the Wgap;default parameters: fasta_swap;optional parameters:fasta_download,fasta_swap,lastdb,last_train,multiz\n')
        self.parser.add_argument('-p', '--parallel', default=1, type=int, help='the number of parallel subprocess')
        self.parser.add_argument('-t', '--threads', type=int, default=1, help='the number of threads for each subprocess')
        self.parser.add_argument('-c', '--configure', help='configure file')
        self.parser.add_argument('-a', '--accession_list', help='accession list', default="")
        self.parser.add_argument('-r', '--ref_fa', help='reference fasta file',required=True)
        self.parser.add_argument('-e', '--exclude_list', default="", help='exclude list')
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


class Wgap_func(object):
    def __init__(self, args):
        # load args
        self.args = args
        self.begin = self.args.begin
        self.parallel = self.args.parallel
        self.threads = self.args.threads
        self.configure = self.args.configure
        self.accession_list = self.args.accession_list
        self.ref_fa = self.args.ref_fa
        self.exclude_list = self.args.exclude_list
        self.out_dir = self.args.out_dir
        self.lastdb = self.args.lastdb
        self.last_train = self.args.last_train
        self.lastal = self.args.lastal
        self.iqtree = self.args.iqtree

        self.acc_list = []
        self.run_dict = {}
        self.result_list = []

        # preprocess
        self.preprocess()
        self.run()

    def preprocess(self):
        # get accession list
        if self.accession_list != "":
            with open(self.accession_list) as f:
                lines = f.readlines()
            for line in lines:
                self.acc_list.append(line.strip())

        # get a fna_dict of one-to-one correspondence between fatsa path and name
        fasta = []
        name = []
        with open(self.configure, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith('#') or line.startswith('\n'):
                continue
            else:
                fasta.append(line.split('\t')[0].strip())
                name.append(line.split('\t')[1].strip())
        fna_dict = dict(zip(name, fasta))

        # get not run list
        not_run_list = list(self.ref_fa.split())
        if self.exclude_list != "":
            with open(self.exclude_list, 'r') as f:
                lines = f.readlines()
            for line in lines:
                not_run_list.append(line.strip())

        # get run list
        for key in fna_dict.keys():
            if key in not_run_list:
                continue
            else:
                self.run_dict[key] = fna_dict[key]

        # mkdir out_dir
        if os.path.isabs(self.out_dir):
            if self.out_dir.endswith('/'):
                self.out_dir = self.out_dir[:-1]
            else:
                pass
        else:
            self.out_dir = os.path.abspath(self.out_dir)
        for i in ["/00_assembly_fasta", "/01_lastdb", "/02_last_train", "/03_lastal", "/04_sort", "/05_multiz", "/06_maf2lst_fa", "/07_iqtree"]:
            sp.Popen("mkdir -p {}{}".format(self.out_dir,i),shell=True).wait()

    def run(self):
        if self.begin == "fasta_download":
            if self.run_fasta_download() == 0:
                print("all fasta download process finished")
                sys.exit(0)
            else:
                print("fasta download process failed")
                sys.exit(1)
        elif self.begin == "lastdb":
            if self.run_lastdb() == 0:
                if self.run_last_train_lastal_sort() == 0:
                    if self.run_multiz() == 0:
                        if self.run_maf2lst_fa() == 0:
                            if self.run_iqtree() == 0:
                                print("all process finished")
                                sys.exit(0)
                            else:
                                print("iqtree process failed")
                                sys.exit(1)
                        else:
                            print("maf2lst_fa process failed")
                            sys.exit(1)
                    else:
                        print("multiz process failed")
                        sys.exit(1)
                else:
                    print("last_train_lastal_sort process failed")
                    sys.exit(1)
            else:
                print("lastdb process failed")
                sys.exit(1)
        elif self.begin == "last_train":
            if self.run_last_train_lastal_sort() == 0:
                if self.run_multiz() == 0:
                    if self.run_maf2lst_fa() == 0:
                        if self.run_iqtree() == 0:
                            print("all process finished")
                            sys.exit(0)
                        else:
                            print("iqtree process failed")
                            sys.exit(1)
                    else:
                        print("maf2lst_fa process failed")
                        sys.exit(1)
                else:
                    print("multiz process failed")
                    sys.exit(1)
            else:
                print("last_train_lastal_sort process failed")
                sys.exit(1)
        elif self.begin == "multiz":
            if self.run_multiz() == 0:
                if self.run_maf2lst_fa() == 0:
                    if self.run_iqtree() == 0:
                        print("all process finished")
                        sys.exit(0)
                    else:
                        print("iqtree process failed")
                        sys.exit(1)
                else:
                    print("maf2lst_fa process failed")
                    sys.exit(1)
            else:
                print("multiz process failed")
                sys.exit(1)
        else:
            if self.run_fasta_swap() == 0:
                if self.run_lastdb() == 0:
                    if self.run_last_train_lastal_sort() == 0:
                        if self.run_multiz() == 0:
                            if self.run_maf2lst_fa() == 0:
                                if self.run_iqtree() == 0:
                                    print("all process finished")
                                    sys.exit(0)
                                else:
                                    print("iqtree process failed")
                                    sys.exit(1)
                            else:
                                print("maf2lst_fa process failed")
                                sys.exit(1)
                        else:
                            print("multiz process failed")
                            sys.exit(1)
                    else:
                        print("last_train_lastal_sort process failed")
                        sys.exit(1)
                else:
                    print("lastdb process failed")
                    sys.exit(1)
            else:
                print("fasta_swap process failed")
                sys.exit(1)


    def run_fasta_download(self):
        # download fasta
        result_list = []
        pool = mp.Pool(self.parallel)
        for i in range(len(self.acc_list)):
            pool.apply_async(self.fasta_download_func, args=(self.acc_list[i],), callback=lambda x: result_list.append(x))
        pool.close()
        pool.join()
        if all([i == 0 for i in result_list]) and len(self.acc_list) == len(result_list):
            return 0
        else:
            return 1

    def run_fasta_swap(self):
        # run fasta_swap
        result_list = []
        pool = mp.Pool(self.parallel)
        for i in range(len(self.run_dict)):
            pool.apply_async(self.fasta_swap_func, args=(list(self.run_dict.keys())[i],), callback=lambda x: result_list.append(x))
        pool.close()
        pool.join()
        if all([i == 0 for i in result_list]) and len(self.run_dict) == len(result_list):
            return 0
        else:
            return 1

    def run_lastdb(self):
        # run lastdb
        if self.lastdb_func() == 0:
            return 0
        else:
            return 1

    def run_last_train_lastal_sort(self):
        # run last_train, lastal, sort
        result_list = []
        pool = mp.Pool(self.parallel)
        for i in range(len(self.run_dict)):
            pool.apply_async(self.last_train_func, args=(list(self.run_dict.keys())[i],), callback=lambda x: result_list.append(x))
        pool.close()
        pool.join()
        if all([i == 0 for i in result_list]) and len(self.run_dict) == len(result_list):
            return 0
        else:
            return 1

    def run_maf2lst_fa(self):
        # run maf2lst_fa
        if self.maf2lst_fa_func() == 0:
            return 0
        else:
            return 1

    def run_multiz(self):
        # run multiz
        result_list = []
        sp.Popen(shlex.split("cd {}{}".format(self.out_dir, "/04_sort/")))
        os.chdir(r"{}{}".format(self.out_dir, "/04_sort/"))

        maf_list = ["{}{}".format(self.out_dir, "/04_sort/")+i+".maf" for i in self.run_dict.keys()]
        # from small to big
        self.maf_list_sort_by_size = sorted(maf_list, key=lambda x: os.path.getsize(x))
        pool = mp.Pool(self.parallel)
        n = 2
        while True:
            if len(result_list) == len(self.run_dict)-1:
                pool.apply_async(self.multiz_func, args=(self.maf_list_sort_by_size[0], self.maf_list_sort_by_size[1], str(n)+'.maf',), callback=lambda x: result_list.append(x))
                break
            else:
                if len(self.maf_list_sort_by_size) <= 1:
                    pass
                else:
                    pool.apply_async(self.multiz_func, args=(self.maf_list_sort_by_size[0], self.maf_list_sort_by_size[1], str(n)+'.maf',), callback=lambda x: result_list.append(x))
                    n += 1
                    self.maf_list_sort_by_size = self.maf_list_sort_by_size[2:]
        pool.close()
        pool.join()

        if all([i == 0 for i in result_list]) and len(self.run_dict)-1 == len(result_list):
            return 0
        else:
            return 1

    def run_iqtree(self):
        # run iqtree
        if self.iqtree_func() == 0:
            return 0
        else:
            return 1

    def fasta_download_func(self, acc):
        # download fasta
        sp.Popen("cd {}".format(self.out_dir),shell=True)
        with open("{}_logfile:".format(acc), 'w', encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directiry:{}\n".format(self.out_dir))
            f.write("\nrunning command line:rsync --copy-links recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/{}* {}\n".format(acc[0:3], acc[4:7], acc[7:10], acc[10:13],acc, self.out_dir))
        f.close()

        cmd_0 = "rsync --copy-links recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/{}* {}".format(acc[0:3], acc[4:7], acc[7:10], acc[10:13],acc, self.out_dir)
        sp_0 = sp.Popen(shlex.split(cmd_0), stdout=sp.PIPE, stderr=sp.PIPE)
        sp_0.communicate()

        if sp_0.returncode == 0:
            with open("{}_logfile:".format(acc), 'a', encoding='utf-8') as f:
                f.write(sp_0.stdout.read().decode('utf-8'))
                f.write("\nfasta_download finished\n")
            f.close()
            return 0
        else:
            with open("error_{}:".format(acc), 'a', encoding='utf-8') as f:
                f.write("\nfasta_download failed\n" + "\nfasta download failed message:\n")
                f.write(sp_0.stderr.read().decode('utf-8'))
            f.close()
            return 1

    def fasta_swap_func(self, name):
        # swap fasta
        fasta = self.run_dict[name]
        sp.Popen(shlex.split("cd {}{}".format(self.out_dir, "/00_assembly_fasta/")))
        os.chdir("{}/00_assembly_fasta/".format(self.out_dir))

        with open("{}_logfile".format(name), "w", encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}{}\n".format(self.out_dir, "/00_assembly_fasta/"))
            f.write(" ".join(['\nrunning command line:zcat', '{}'.format(fasta), '|', 'awk', '\'{print $1}\'', '|', 'sed', '\'s/>/>{}_/g\''.format(name), '>', '{}.fa\n'.format(name)]))
        f.close()

        sp_0 = sp.Popen(["zcat", "{}".format(fasta)], stdout=sp.PIPE)
        sp_1 = sp.Popen(shlex.split("awk '{print $1}'"), stdin=sp_0.stdout, stdout=sp.PIPE)
        sp_0.stdout.close()
        sp_2 = sp.Popen(shlex.split("sed 's/>/>{}_/g'".format(name)), stdin=sp_1.stdout, stdout=sp.PIPE)
        sp_1.stdout.close()

        with open('{}.fa'.format(name), 'w', encoding='utf-8') as f:
            for line in iter(sp_2.stdout.readline, b''):
                f.write(line.decode('utf-8'))
        f.close()
        sp_2.communicate()
        sp_2.stdout.close()

        if sp_2.returncode == 0:
            with open("{}_logfile".format(name), "a", encoding='utf-8') as f:
                f.write("\nfasta_swap successfully finished\n")
                f.write("\nrunning command line:samtools faidx {}.fa\n".format(name))
            f.close()

        sp_3 = sp.Popen(shlex.split("samtools faidx {}.fa".format(name)), stdout=sp.PIPE, stderr=sp.PIPE)
        with open("{}.fa.fai".format(name), "w", encoding='utf-8') as f:
            for line in iter(sp_3.stdout.readline, b''):
                f.write(line.decode('utf-8'))
        sp_3.communicate()
        sp_3.stdout.close()

        if sp_3.returncode == 0:
            with open("{}_logfile".format(name), "a", encoding='utf-8') as f:
                f.write("\nsamtools faidx successfully finished\n")
            f.close()
            return 0
        else:
            with open("error_{}".format(name), "w", encoding='utf-8') as f:
                f.write("\nsamtools faidx failed\n" + "\nsamtools faidx failed message:\n")
                f.write(sp_3.stderr.read().decode('UTF-8').strip() + '\n')
            f.close()
            return 1

    def lastdb_func(self):
        sp.Popen(shlex.split("cd {}{}".format(self.out_dir, "/01_lastdb/")))
        os.chdir(r"{}/01_lastdb/".format(self.out_dir))

        with open("logfile", "w", encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}{}\n".format(self.out_dir, "/01_lastdb/"))
            f.write("\nrunning command line:lastdb -P {} ".format(self.parallel * self.threads) + self.lastdb + " " + self.ref_fa + '_db ' + self.out_dir + '/00_assembly_fasta/' + self.ref_fa + '.fa\n')
        f.close()

        cmd_0 = "lastdb -P {} ".format(self.parallel * self.threads) + self.lastdb + " " + self.ref_fa + '_db ' + self.out_dir + '/00_assembly_fasta/' + self.ref_fa + '.fa'
        sp_0 = sp.Popen(cmd_0, stdout=sp.PIPE, stderr=sp.PIPE)
        sp_0.communicate()
        if sp_0.returncode == 0:
            with open("logfile".format(self.out_dir), "a", encoding='utf-8') as f:
                f.write("\nlastdb successfully finished\n")
            f.close()
            return 0
        else:
            sp_1 = sp.Popen(shlex.split("lastdb --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            with open("errorfile", "w", encoding='utf-8') as f:
                f.write("\nlastdb failed\n" + "lastdb failed message:\n")
                f.write(sp_0.stderr.read().decode('UTF-8'))
                f.write("\nlastdb help:\n")
                f.write(sp_1.stdout.read().decode('UTF-8'))
                f.write(sp_1.stderr.read().decode('UTF-8'))
            f.close()
            return 1

    def last_train_func(self, name):
        sp.Popen(shlex.split("cd {}{}".format(self.out_dir,"/02_last_train/")))
        os.chdir(r"{}{}".format(self.out_dir, "/02_last_train/"))

        with open("{}_logfile".format(name), "w", encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}{}\n".format(self.out_dir, "/02_last_train/"))
            f.write("\nrunning command line:last-train -P {} {} {}/01_lastdb/{}_db {}/00_assembly_fasta/{}.fa > {}/02_last_train/{}.mat\n".format(self.threads, self.last_train, self.out_dir, self.ref_fa, self.out_dir, name, self.out_dir, name))
        f.close()

        cmd_0 = "last-train -P {} ".format(self.threads) + self.last_train + " " + self.out_dir + "/01_lastdb/" + self.ref_fa + "_db " + self.out_dir + "/00_assembly_fasta/" + name + ".fa"

        with open("{}.mat".format(name), "w") as f:
            sys.stdout = f
            sp_0 = sp.Popen(shlex.split(cmd_0), stdout=sys.stdout, stderr=sp.PIPE)
            sp_0.communicate()
        sys.stdout.close()
        sp_0.stdout.close()
        f.close()

        # check last-train return code
        if sp_0.returncode == 0:
            with open("{}_logfile".format(name), "a", encoding='utf-8') as f:
                f.write("\nlast-train successfully finished\n")
            f.close()
            return self.lastal_func(name)
        else:
            with open("error_{}".format(name), "w",encoding='utf-8') as f:
                f.write("\nlast-train failed\n" + "last-train error message:\n")
                f.write(sp_0.stderr.read().decode('UTF-8'))
            f.close()
            sp_0.stderr.close()
            return 1

    def lastal_func(self, name):
        sp.Popen(shlex.split("cd {}{}".format(self.out_dir, "/03_lastal/")))
        os.chdir(r"{}{}".format(self.out_dir, "/03_lastal/"))

        with open("{}_logfile".format(name), "w") as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}{}\n".format(self.out_dir, "/03_lastal/"))
            f.write("\nrunning command line:lastal -P {} {} {}/02_last_train/{}.mat {}/01_lastdb/{}_db {}/00_assembly_fasta/{}.fa | last-split -f MAF+ > {}/03_lastal/{}.maf\n".format(self.threads, self.lastal, self.out_dir, name, self.out_dir, self.ref_fa, self.out_dir, name, self.out_dir, name))
        f.close()

        cmd_0 = "lastal -P {} ".format(self.threads) + self.lastal + " " + self.out_dir + "/02_last_train/" + name + ".mat " + self.out_dir + "/01_lastdb/" + self.ref_fa + "_db " + self.out_dir + "/00_assembly_fasta/" + name + ".fa "
        cmd_1 = "last-split " + "-f " + "MAF+"


        with open("{}.maf".format(name), "w", encoding='utf-8') as f:
            sys.stdout = f
            sp_0 = sp.Popen(shlex.split(cmd_0), stdout=sp.PIPE, stderr=sp.PIPE)
            sp_1 = sp.Popen(shlex.split(cmd_1), stdin=sp_0.stdout, stdout=sys.stdout, stderr=sp.PIPE)
            sp_1.communicate()
        sys.stdout.close()
        sp_1.stdout.close()
        sp_0.stdout.close()
        f.close()

        # check sp_1 return code
        if sp_1.returncode == 0:
            with open("{}_logfile".format(name), "a", encoding='utf-8') as f:
                f.write("\nlastal and last-split successfully finished\n")
            f.close()
            return self.sort_func(name)
        else:
            with open("error_{}".format(name), "w", encoding='utf-8') as f:
                f.write("\nerror:lastal and last-split failed\n" + "lastal error message:\n")
                f.write(sp_0.stderr.read().decode('UTF-8'))
                f.write("\nlast-split error message:\n")
                f.write(sp_1.stderr.read().decode('UTF-8'))
            f.close()
            sp_0.stderr.close()
            sp_1.stderr.close()
            return 1

    def sort_func(self, name):
        sp.Popen(shlex.split("cd {}{}".format(self.out_dir, "/04_sort/")))
        os.chdir(r"{}{}".format(self.out_dir, "/04_sort/"))

        with open("{}_logfile".format(name), "w") as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}{}\n".format(self.out_dir, "/04_sort/"))
            f.write("\nrunning command line:maf-swap {}/03_lastal/{}.maf |last-split |maf-swap |maf-sort|grep -v '#'|sed '1i\#maf version=1.0 scoring=last' > {}/04_sort/{}.maf\n".format(self.out_dir, name, self.out_dir, name))
        f.close()

        cmd_0 = "maf-swap {}/03_lastal/{}.maf".format(self.out_dir, name)
        cmd_1 = "last-split"
        cmd_2 = "maf-swap"
        cmd_3 = "maf-sort"

        with open("{}.maf".format(name), "w", encoding='utf-8') as f:
            sys.stdout = f
            sys.stdout.write("##maf version=1 scoring=last\n")
            sp_0 = sp.Popen(shlex.split(cmd_0), stdout=sp.PIPE, stderr=sp.PIPE)
            sp_1 = sp.Popen(shlex.split(cmd_1), stdin=sp_0.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
            sp_2 = sp.Popen(shlex.split(cmd_2), stdin=sp_1.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
            sp_3 = sp.Popen(shlex.split(cmd_3), stdin=sp_2.stdout, stdout=sys.stdout, stderr=sp.PIPE)
            sp_3.communicate()
        sys.stdout.close()
        sp_3.stdout.close()
        sp_2.stdout.close()
        sp_1.stdout.close()
        sp_0.stdout.close()
        f.close()

        # check sp_3 return code
        if sp_3.returncode == 0:
            with open("{}_logfile".format(name), "a", encoding='utf-8') as f:
                f.write("\nmaf-sort successfully finished\n")
            f.close()
            return 0
        else:
            with open("error_{}".format(name), "w", encoding='utf-8') as f:
                f.write("\nerror:maf-sort failed\n" + "maf-swap error message:\n")
                f.write(sp_0.stderr.read().decode('UTF-8'))
                f.write("\nlast-split error message:\n")
                f.write(sp_1.stderr.read().decode('UTF-8'))
                f.write("\nmaf-swap error message:\n")
                f.write(sp_2.stderr.read().decode('UTF-8'))
                f.write("\nmaf-sort error message:\n")
                f.write(sp_3.stderr.read().decode('UTF-8'))
            f.close()
            sp_0.stderr.close()
            sp_1.stderr.close()
            sp_2.stderr.close()
            sp_3.stderr.close()
            return 1

    def multiz_func(self,left,right,name):
        sp.Popen(shlex.split("cd {}{}".format(self.out_dir, "/05_multiz/")))
        os.chdir(r"{}{}".format(self.out_dir, "/05_multiz/"))

        with open("{}.logfile".format(name), "w") as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}{}\n".format(self.out_dir, "/05_multiz/"))
            f.write("\nrunning command line:multiz {} {} 0 {}_1 {}_2 > {}\n".format(left,right,name,name,name))
        f.close()

        cmd_0 = "multiz {} {} 0 {}_1 {}_2".format(left,right,name,name)
        with open (name,'w',encoding='utf-8') as f:
            sys.stdout = f
            sp_0 = sp.Popen(shlex.split(cmd_0), stdout=sys.stdout, stderr=sp.PIPE)
            sp_0.communicate()
        sys.stdout.close()
        sp_0.stdout.close()
        f.close()

        # check sp_0 return code
        if sp_0.returncode == 0:
            with open("{}.logfile".format(name), "a", encoding='utf-8') as f:
                f.write("\nmultiz successfully finished\n")
            f.close()
            self.maf_list_sort_by_size.append(name)
            return 0
        else:
            with open("error_{}".format(name), "w", encoding='utf-8') as f:
                f.write("\nerror:multiz failed\n" + "multiz error message:\n")
                f.write(sp_0.stderr.read().decode('UTF-8'))
            f.close()
            sp_0.stderr.close()
            return 1

    def maf2lst_fa_func(self):
        maf_file = "{}/05_multiz/{}.maf".format(self.out_dir, len(self.run_dict))
        lst_file = "{}/06_maf2lst_fa/{}_{}.lst".format(self.out_dir, self.ref_fa, len(self.run_dict))
        fasta_file = "{}/06_maf2lst_fa/{}_{}.fa".format(self.out_dir, self.ref_fa, len(self.run_dict))

        sp.Popen(shlex.split("cd {}{}".format(self.out_dir, "/06_maf2lst_fa/")))
        os.chdir(r"{}{}".format(self.out_dir, "/06_maf2lst_fa/"))

        cmd_0 = 'csplit %s /^a/ -n8 -s {*} -f block_' %(maf_file)
        sp_0 = sp.Popen(shlex.split(cmd_0))
        sp_0.wait()

        for maf_block_file in [f for f in os.listdir('.') if f.startswith('block')]:
            with open(maf_block_file, 'r', encoding='utf-8') as maf_block:
                maf_block_startswith_s = [line for line in maf_block.readlines() if line.startswith('s')]

            if len(maf_block_startswith_s) == len(self.run_dict) + 1:
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

    def iqtree_func(self):
        fasta_file = "{}/06_maf2lst_fa/{}_{}.fa".format(self.out_dir, self.ref_fa, len(self.run_dict))
        sp.Popen(shlex.split("cd {}{}".format(self.out_dir, "/07_iqtree/")))
        os.chdir(r"{}{}".format(self.out_dir, "/07_iqtree/"))

        with open("iqtree_logfile", "w", encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write("\nrunning command line:iqtree -s {} -nt {} {}\n".format(fasta_file, self.threads, self.iqtree))
        f.close()

        cmd_0 = "iqtree -s {} -nt {} {}".format(fasta_file, self.threads, self.iqtree)
        sp_0 = sp.Popen(shlex.split(cmd_0), stdout=sp.PIPE, stderr=sp.PIPE)
        sp_0.wait()

        if sp_0.returncode == 0:
            with open("iqtree_logfile", "a", encoding='utf-8') as f:
                f.write("\n{} iqtree successfully finished\n".format(len(self.run_dict)))
                f.write(sp_0.stdout.read().decode('UTF-8'))
            f.close()
            return 0
        else:
            with open("iqtree_error_log", "w", encoding='utf-8') as f:
                f.write("\nerror:iqtree failed\n" + "iqtree error message:\n")
                f.write(sp_0.stderr.read().decode('UTF-8'))
            f.close()
            sp_0.stderr.close()
            return 1



if __name__ == '__main__':
    args = Wgap_help().args
    Wgap_func(args)

