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
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


_version = "1.0.0"

class wga_in_one_step(object):
    def __init__(self):
        description = 'wga_in_one_step{}: whole genome alignment in one step with last'.format(_version)
        self.parser = argparse.ArgumentParser(
            prog='wga_in_one_step',
            description=description)
        self.parser.add_argument('-t', '--threads', type=int, default=1, help='number of threads')
        self.parser.add_argument('-c', '--configure', help='configure file path', required=False)
        self.parser.add_argument('-f', '--fa_dir', help='fasta file directory',default=os.getcwd())
        self.parser.add_argument('-r', '--ref_fa', help='reference fasta file')
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
        self.ref_fa = self.args.ref_fa
        self.exclude_fa = self.args.exclude_fa
        self.out_dir = self.args.out_dir
        self.lastdb = self.args.lastdb
        self.last_train = self.args.last_train
        self.lastal = self.args.lastal

        self.fasta_swap()
        # self.lastdb_func()
        # self.last_train_func()
        # self.lastal_func()

    def fasta_swap(self):
        sp.Popen(shlex.split("mkdir -p {}/00_assembly_fasta".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/00_assembly_fasta/".format(self.out_dir)))
        os.chdir(r"{}00_assembly_fasta/".format(self.out_dir))
        with open(self.configure, 'r') as f:
            lines = f.readlines()
        for line in lines:
            self.fna_gz = line.split('\t')[0].strip()
            self.name = line.split('\t')[1].strip()

            sys.stdout = open("{}/00_assembly_fasta/{}_logfile".format(self.out_dir,self.name), "w")
            sys.stdout.write("\nargparse:{}\n".format(self.args))
            sys.stdout.write("\nrunning directory:{}\n".format(self.out_dir))
            sys.stdout.write(" ".join(['\nrunning command line:', 'cat', '"{}"'.format(self.fna_gz), '|', 'gunzip', '|', 'awk', '-F', '\'\t\'', '\'{print $1}\'', '|', 'sed', '\'s/>/>{}_/g\''.format(self.name), '>', '{}/00_assembly_fasta/{}.fa\n'.format(self.out_dir,self.name)]))
            sys.stdout.flush()

            self.sp_ungzip = sp.Popen(['cat', '"{}"'.format(self.fna_gz), '|', 'gunzip', '|', 'awk', '-F', '\t', '{print $1}', '|', 'sed', 's/>/>{}_/g'.format(self.name)], stdout=sp.PIPE)
            for fa_line in self.sp_ungzip.stdout:
                fa_line = line.decode('utf-8')
                with open('{}/00_assembly_fasta/{}.fa'.format(self.out_dir,self.name), 'a') as fa:
                    fa.write(line)
            self.sp_ungzip.stdout.flush()
            self.sp_ungzip.communicate()
            self.sp_ungzip.stdout.close()

            if self.sp_ungzip.returncode == 0:
                sys.stdout.write("\nfasta_swap_1 successfully finished\n")
                sys.stdout.write("\nrunning command line:samtools faidx -@ {} {}/00_assembly_fasta/{}.fa\n".format(self.threads, self.out_dir, self.name))
                sys.stdout.flush()
                self.sp_faidx = sp.Popen(shlex.split("samtools faidx {}/00_assembly_fasta/{}.fa".format(self.out_dir,self.name)), stdout=sp.PIPE)
                for fai_line in self.sp_faidx.stdout:
                    fai_line = line.decode('utf-8')
                    with open('{}/00_assembly_fasta/{}.fa.fai'.format(self.out_dir,self.name), 'a') as fai:
                        fai.write(line)
                self.sp_faidx.stdout.flush()
                self.sp_faidx.communicate()
                self.sp_faidx.stdout.close()


                if self.sp_faidx.returncode == 0:
                    sys.stdout.write("\nfasta_swap_2 successfully finished\n")
                else:
                    sys.stdout.write("\nfasta_swap_2 failed\n")
            else:
                sys.stdout.write("\nfasta_swap_1 failed\n")
            sys.stdout.flush()
            sys.stdout.close()


    def lastdb_func(self):
        sp.Popen(shlex.split("mkdir -p {}/00_lastdb".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/00_lastdb/".format(self.out_dir)))

        sys.stdout = open("{}/00_lastdb/logfile".format(self.out_dir), "w")
        sys.stdout.write("\nargparse:{}\n".format(self.args))
        sys.stdout.write("\nrunning directory:{}\n".format(self.out_dir))
        sys.stdout.write("\nrunning command line:lastdb -P {} {}\n".format(self.threads, self.lastdb))
        sys.stdout.flush()

        os.chdir(r"{}/00_lastdb/".format(self.out_dir))
        self.sp = sp.Popen(shlex.split("lastdb -P {} {}".format(self.threads, self.lastdb)), stdout=sp.PIPE, stderr=sp.PIPE)
        # wiat for the process to finish
        self.sp.communicate()
        # check lastdb return code
        if self.sp.returncode == 0:
            sys.stdout = open("{}/00_lastdb/logfile".format(self.out_dir), "a")
            sys.stdout.write("\nlastdb successfully finished\n")
            sys.stdout.flush()
        else:
            sys.stdout = open("{}/00_lastdb/err.file".format(self.out_dir), "a")
            sys.stdout.write("\nerror:lastdb failed\n"+"lastdb error message:\n")
            for stdout_line in iter(self.sp.stderr.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            sys.stdout.flush()

            self.sp_help = sp.Popen(shlex.split("lastdb --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            sys.stdout = open("{}/00_lastdb/err.file".format(self.out_dir), "a")
            for stdout_line in iter(self.sp_help.stderr.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            sys.stdout.flush()



    def last_train_func(self):
        sp.Popen(shlex.split("mkdir -p {}/01_last_train".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/01_last_train/".format(self.out_dir)))

        sys.stdout = open("{}/01_last_train/logfile".format(self.out_dir), "w")
        sys.stdout.write("\nargparse:{}\n".format(self.args))
        sys.stdout.write("\nrunning directory:{}\n".format(self.out_dir))
        sys.stdout.write("\nrunning command line:last-train {} {}\n".format(self.last_train, self.ref_fa))
        sys.stdout.flush()

        os.chdir(r"{}/01_last_train/".format(self.out_dir))
        self.sp = sp.Popen(shlex.split("last-train {} {}".format(self.last_train, self.ref_fa)), stdout=sp.PIPE, stderr=sp.PIPE)
        # wiat for the process to finish
        self.sp.communicate()
        # check lastdb return code
        if self.sp.returncode == 0:
            sys.stdout = open("{}/01_last_train/logfile".format(self.out_dir), "a")
            sys.stdout.write("\nlast-train successfully finished\n")
            sys.stdout.flush()
        else:
            sys.stdout = open("{}/01_last_train/logfile".format(self.out_dir), "a")
            sys.stdout.write("\nerror:last-train failed\n"+"last-train error message:\n")
            for stdout_line in iter(self.sp.stderr.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            sys.stdout.flush()


    def lastal_func(self):
        sp.Popen(shlex.split("mkdir -p {}/02_lastal".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/02_lastal/".format(self.out_dir)))

        sys.stdout = open("{}/02_lastal/logfile".format(self.out_dir), "w")
        sys.stdout.write("\nargparse:{}\n".format(self.args))
        sys.stdout.write("\nrunning directory:{}\n".format(self.out_dir))
        sys.stdout.write("\nrunning command line:lastal -P {} {} {}\n".format(self.threads, self.lastal, self.ref_fa, self.query_fa))
        sys.stdout.flush()

        os.chdir(r"{}/02_lastal/".format(self.out_dir))
        self.sp = sp.Popen(shlex.split("lastal -P {} {} {} {}".format(self.threads, self.lastal, self.ref_fa, self.query_fa)), stdout=sp.PIPE, stderr=sp.PIPE)
        # wiat for the process to finish
        self.sp.communicate()
        # check lastdb return code
        if self.sp.returncode == 0:
            sys.stdout = open("{}/02_lastal/logfile".format(self.out_dir), "a")
            sys.stdout.write("\nlastal successfully finished\n")
            sys.stdout.flush()
        else:
            sys.stdout = open("{}/02_lastal/err.file".format(self.out_dir), "a")
            sys.stdout.write("\nerror:lastal failed\n"+"lastal error message:\n")
            for stdout_line in iter(self.sp.stderr.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            sys.stdout.flush()



if __name__ == '__main__':
    wga_in_one_step()
