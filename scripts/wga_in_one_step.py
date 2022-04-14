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


_version = "1.0.0"

class wga_in_one_step(object):
    def __init__(self):
        description = 'wga_in_one_step: whole genome alignment in one step with last\n'\
                      'https://github.com/hyBio/L.guttatus/blob/master/scripts/wga_in_one_step.py\n'
        self.parser = argparse.ArgumentParser(
            prog='wga_in_one_step',
            description=description)
        self.parser.add_argument('-b', '--begin', default="fasta_swap", help='the begin subprocess of the pipeline\n'
                                                                             'default: fasta_swap\n'
                                                                             'optional parameters:fasta_download\tfasta_swap\tlastdb\tlast_train\tlastal\tsort\tmultiz',required=True)
        self.parser.add_argument('-p', '--parallel', default=1, type=int, help='the number of parallel subprocess')
        self.parser.add_argument('-t', '--threads', type=int, default=1, help='the number of threads for each subprocess')
        self.parser.add_argument('-c', '--configure', help='configure file path')
        self.parser.add_argument('-a', '--accession_list', help='accession number list file path', default="")
        self.parser.add_argument('-f', '--fa_dir', help='fasta file directory',default=os.getcwd())
        self.parser.add_argument('-r', '--ref_fa', help='reference fasta file')
        self.parser.add_argument('-e', '--exclude_fa', default="", help='exclude fasta file name')
        self.parser.add_argument('-o', '--out_dir', help='output directory', default=os.getcwd())
        self.parser.add_argument('-s1', '--lastdb', help='lastdb argparse',type=str,default='')
        self.parser.add_argument('-s2', '--last_train', help='last_train argparse',type=str,default='')
        self.parser.add_argument('-s3', '--lastal', help='lastal argparse',type=str,default='')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + _version)

        if len(sys.argv) == 1:
            self.parser.print_help()
            sys.exit(1)
        else:
            self.args = self.parser.parse_args()
            self.whole_genome_alignment = whole_genome_alignment(self.args)




class whole_genome_alignment(object):

    def __init__(self,args):
        self.args = args
        self.parallel = self.args.parallel
        self.begin = self.args.begin
        self.threads = self.args.threads
        self.configure = self.args.configure
        self.fa_dir = self.args.fa_dir
        self.ref_fa = self.args.ref_fa
        self.exclude_fa = self.args.exclude_fa
        self.out_dir = self.args.out_dir
        self.lastdb = self.args.lastdb
        self.last_train = self.args.last_train
        self.lastal = self.args.lastal

        self.preprocessing()
        self.running()

    def preprocessing(self):
        self.accession_list = []
        self.fna_gz_list = []
        self.name_list = []
        self.exclude_fa_list = []

        if self.accession_list != "":
            with open(self.accession_list) as f:
                lines = f.readlines()
            for line in lines:
                self.accession_list.append(line.strip())

        if self.exclude_fa != "":
            with open(self.exclude_fa, 'r') as f:
                lines = f.readlines()
            for line in lines:
                self.exclude_fa_list.append(line.strip())

        with open(self.configure, 'r') as f:
            lines = f.readlines()
        for line in lines:
            self.fna_gz_list.append(line.split('\t')[0].strip())
            self.name_list.append(line.split('\t')[1].strip())
        self.fasta_swap_name_list = self.name_list
        self.name_list.remove(self.ref_fa)
        for i in self.exclude_fa_list:
            self.name_list.remove(i)
        self.last_train_name_list = self.name_list
        self.lastal_name_list = self.name_list
        self.sort_name_list = self.name_list


    def running(self):
        if self.begin == "fasta_download":
            self.fasta_download_run()
        elif self.begin == "fasta_swap":
            self.fasta_swap_run()
            self.lastdb_run()
            self.last_train_run()
            self.lastal_run()
            self.sort_run()
        elif self.begin == "lastdb":
            self.lastdb_run()
            self.last_train_run()
            self.lastal_run()
            self.sort_run()
        elif self.begin == "last_train":
            self.last_train_run()
            self.lastal_run()
            self.sort_run()
        elif self.begin == "lastal":
            self.lastal_run()
            self.sort_run()
        elif self.begin == "sort":
            self.sort_run()
        else:
            print("Incorrect parameter '{}' given to whole_genome_alignmen, please check the help".format(self.begin))
            sys.exit(1)


    def fasta_download_func(self,accession):
        self.accession = accession
        with open("{}_logfile:".format(self.accession),'w',encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directiry:{}\n".format(self.fa_dir))
            f.write("\nrsync --copy-links recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/{}* {}\n".format(self.accession[0:3],self.accession[4:7],self.accession[7:10],self.accession[10:13],self.accession,self.fa_dir,self.accession))
        f.close()

        self.sp_fasta_download = sp.Popen(shlex.split("rsync --copy-links recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/{}* {}".format(self.accession[0:3],self.accession[4:7],self.accession[7:10],self.accession[10:13],self.accession,self.fa_dir)),stdout=sp.PIPE,stderr=sp.PIPE)
        self.sp_fasta_download.communicate()

        if self.sp_fasta_download.returncode == 0:
            with open("{}_logfile:".format(self.accession),'a',encoding='utf-8') as f:
                f.write(self.sp_fasta_download.stdout.read().decode('utf-8'))
                f.write("\nfasta_download finished\n")
            f.close()
            return(0)
        else:
            with open("error_{}:".format(self.accession),'a',encoding='utf-8') as f:
                f.write("\nfasta_download failed\n" + "\nfasta download failed message:\n")
                f.write(self.sp_fasta_download.stderr.read().decode('utf-8'))
            f.close()
            return(1)


    def fasta_swap_func(self,fna_gz,fasta_swap_name):
        self.fna_gz = fna_gz
        self.fasta_swap_name = fasta_swap_name

        sp.Popen(shlex.split("mkdir -p {}/00_assembly_fasta".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/00_assembly_fasta/".format(self.out_dir)))
        os.chdir(r"{}00_assembly_fasta/".format(self.out_dir))

        with open("{}/00_assembly_fasta/{}_logfile".format(self.out_dir,fasta_swap_name.name), "w", encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write(" ".join(['\nrunning command line:', 'zcat', '{}'.format(self.fna_gz),'|', 'awk', '\'{print $1}\'', '|', 'sed', '\'s/>/>{}_/g\''.format(self.fasta_swap_name), '>', '{}/00_assembly_fasta/{}.fa\n'.format(self.out_dir,self.fasta_swap_name)]))
        f.close()

        self.sp_ungzip = sp.Popen(["zcat", "{}".format(self.fna_gz)], stdout=sp.PIPE)
        self.sp_awk = sp.Popen(shlex.split("awk '{print $1}'"), stdin=self.sp_ungzip.stdout, stdout=sp.PIPE)
        self.sp_ungzip.stdout.close()
        self.sp_sed = sp.Popen(shlex.split("sed 's/>/>{}_/g'".format(self.fasta_swap_name)), stdin=self.sp_awk.stdout, stdout=sp.PIPE)
        self.sp_awk.stdout.close()
        fasta_swap_1_output = self.sp_sed.communicate()[0]
        with open('{}/00_assembly_fasta/{}.fa'.format(self.out_dir,self.fasta_swap_name), 'a', encoding='utf-8') as fa:
            fa.write(fasta_swap_1_output.decode('UTF-8').strip() + '\n')
        self.sp_sed.stdout.close()

        if self.sp_sed.returncode == 0:
            with open("{}/00_assembly_fasta/{}_logfile".format(self.out_dir,self.fasta_swap_name), "a", encoding='utf-8') as f:
                f.write("\nfasta_swap successfully finished\n")
                f.write("\nrunning command line:samtools faidx {}/00_assembly_fasta/{}.fa\n".format(self.out_dir, self.fasta_swap_name))
            f.close()
            self.sp_faidx = sp.Popen(shlex.split("samtools faidx {}/00_assembly_fasta/{}.fa".format(self.out_dir,self.fasta_swap_name)), stdout=sp.PIPE)
            fasta_swap_2_output = self.sp_faidx.communicate()[0]
            with open('{}/00_assembly_fasta/{}.fa.fai'.format(self.out_dir,self.fasta_swap_name), 'a', encoding='utf-8') as f:
                f.write(fasta_swap_2_output.decode('UTF-8').strip() + '\n')
            self.sp_faidx.stdout.close()

            if self.sp_faidx.returncode == 0:
                with open("{}/00_assembly_fasta/{}_logfile".format(self.out_dir,self.fasta_swap_name), "a", encoding='utf-8') as f:
                    f.write("\nsamtools faidx successfully finished\n")
                f.close()
                return (0)
            else:
                with open("{}/00_assembly_fasta/error_{}".format(self.out_dir,self.fasta_swap_name), "a", encoding='utf-8') as f:
                    f.write("\nsamtools faidx failed\n"+"samtools faidx failed message:\n")
                    f.write(self.sp_faidx.stderr.read().decode('UTF-8'))
                f.close()
                return(1)
        else:
            with open("{}/00_assembly_fasta/error_{}".format(self.out_dir,self.fasta_swap_name), "a", encoding='utf-8') as f:
                f.write("\nsed failed\n"+"sed failed message:\n")
                f.write(self.sp_sed.stderr.read().decode('UTF-8'))
            f.close()
            return(1)


    def lastdb_func(self):
        sp.Popen(shlex.split("mkdir -p {}/01_lastdb".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/01_lastdb/".format(self.out_dir)))

        with open("{}/01_lastdb/logfile".format(self.out_dir), "w", encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write("\nrunning command line:lastdb -P {} ".format(self.parallel*self.threads) + self.lastdb + " " + self.ref_fa + '_db ' + self.out_dir + '/00_assembly_fasta/' + self.ref_fa + '.fa\n')
        f.close()

        os.chdir(r"{}/01_lastdb/".format(self.out_dir))
        sp_lastdb_cmd = "lastdb -P {} ".format(self.parallel*self.threads) + self.lastdb + " " + self.ref_fa + '_db ' + self.out_dir + '/00_assembly_fasta/' + self.ref_fa + '.fa'
        self.sp_lastdb = sp.Popen(shlex.split(sp_lastdb_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
        # wiat for the process to finish
        self.sp_lastdb.communicate()
        # check lastdb return code
        if self.sp_lastdb.returncode == 0:
            with open("{}/01_lastdb/logfile".format(self.out_dir), "a", encoding='utf-8') as f:
                f.write("\nlastdb successfully finished\n")
            f.close()
            return(0)
        else:
            self.sp_lastdb_help = sp.Popen(shlex.split("lastdb --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            with open("{}/01_lastdb/logfile".format(self.out_dir), "a", encoding='utf-8') as f:
                f.write("\nlastdb failed\n"+"lastdb failed message:\n")
                f.write(self.sp_lastdb.stderr.read().decode('UTF-8'))
                f.write("\nlastdb help:\n")
                f.write(self.sp_lastdb_help.stdout.read().decode('UTF-8'))
                f.write(self.sp_lastdb_help.stderr.read().decode('UTF-8'))
            f.close()
            return(1)


    def last_train_func(self,last_train_name):
        self.last_train_name = last_train_name
        sp.Popen(shlex.split("mkdir -p {}/02_last_train".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/02_last_train/".format(self.out_dir)))
        os.chdir(r"{}/02_last_train/".format(self.out_dir))

        with open("{}/02_last_train/{}_logfile".format(self.out_dir,self.last_train_name), "w", encoding='utf-8') as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write("\nrunning command line:last-train -P {} {} {}/01_lastdb/{}_db {}/00_assembly_fasta/{}.fa > {}/02_last_train/{}.mat\n".format(self.threads, self.last_train, self.out_dir, self.ref_fa, self.out_dir, self.last_train_name, self.out_dir, self.last_train_name))
        f.close()

        sp_last_train_cmd = "last-train -P {} ".format(self.threads) + self.last_train + " " + self.out_dir + "/01_lastdb/" + self.ref_fa + "_db " + self.out_dir + "/00_assembly_fasta/" + self.last_train_name + ".fa"
        self.sp_last_train = sp.Popen(shlex.split(sp_last_train_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
        # wiat for the process to finish
        last_train_output = self.sp_last_train.communicate()[0]
        with open("{}/02_last_train/{}.mat".format(self.out_dir, self.last_train_name), "a") as f:
            f.write(last_train_output.decode('UTF-8').strip()+"\n")
        self.sp_last_train.stdout.close()

        # check lastdb return code
        if self.sp_last_train.returncode == 0:
            with open("{}/02_last_train/{}_logfile".format(self.out_dir,self.last_train_name), "a", encoding='utf-8') as f:
                f.write("\nlast-train successfully finished\n")
            f.close()
            return(0)
        else:
            with open("{}/02_last_train/error_{}".format(self.out_dir,self.last_train_name), "a", encoding='utf-8') as f:
                f.write("\nlast-train failed\n" + "last-train error message:\n")
                f.write(self.sp_last_train.stderr.read().decode('UTF-8'))
            f.close()
            self.sp_last_train.stderr.close()
            return(1)


    def lastal_func(self,lastal_name):
        self.lastal_name = lastal_name
        sp.Popen(shlex.split("mkdir -p {}/03_lastal".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/03_lastal/".format(self.out_dir)))
        os.chdir(r"{}/03_lastal/".format(self.out_dir))

        with open("{}/03_lastal/{}_logfile".format(self.out_dir, self.lastal_name), "w") as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write("\nrunning command line:lastal -P {} {} {}/02_last_train/{}.mat {}/01_lastdb/{}_db {}/00_assembly_fasta/{}.fa | last-split -f MAF+ > {}/03_lastal/{}.maf\n".format(self.threads, self.lastal, self.out_dir, self.lastal_name, self.out_dir, self.ref_fa, self.out_dir, self.lastal_name, self.out_dir, self.lastal_name))
        f.close()

        sp_lastal_1_cmd = "lastal -P {} ".format(self.threads) + self.lastal + " " + self.out_dir + "/02_last_train/" + self.lastal_name + ".mat " + self.out_dir + "/01_lastdb/" + self.ref_fa + "_db " + self.out_dir + "/00_assembly_fasta/" + self.lastal_name + ".fa "
        sp_lastal_2_cmd = "last-split " + "-f " + "MAF+"

        self.sp_lastal_1 = sp.Popen(shlex.split(sp_lastal_1_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
        self.sp_lastal_2 = sp.Popen(shlex.split(sp_lastal_2_cmd), stdin=self.sp_lastal_1.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
        sp_lastal_2_output = self.sp_lastal_2.communicate()[0]
        with open("{}/03_lastal/{}.maf".format(self.out_dir, self.lastal_name), "a", encoding='utf-8') as f:
            f.write(sp_lastal_2_output.decode('UTF-8').strip()+"\n")
        f.close()
        self.sp_lastal_1.stdout.close()

        # check lastal_2 return code
        if self.sp_lastal_2.returncode == 0:
            with open("{}/03_lastal/{}_logfile".format(self.out_dir, self.lastal_name), "a", encoding='utf-8') as f:
                f.write("\nlastal and last-split successfully finished\n")
            f.close()
            return(0)
        else:
            with open("{}/03_lastal/error_{}".format(self.out_dir, self.lastal_name), "a", encoding='utf-8') as f:
                f.write("\nerror:lastal and last-split failed\n" + "lastal error message:\n")
                f.write(self.sp_lastal_1.stderr.read().decode('UTF-8'))
                f.write("\nlast-split error message:\n")
                f.write(self.sp_lastal_2.stderr.read().decode('UTF-8'))
            f.close()
            self.sp_lastal_1.stderr.close()
            self.sp_lastal_2.stderr.close()
            return(1)


    def sort_func(self,sort_name):
        self.sort_name = sort_name
        sp.Popen(shlex.split("mkdir -p {}/04_sort".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/04_sort/".format(self.out_dir)))
        os.chdir(r"{}/04_sort/".format(self.out_dir))

        with open("{}/04_sort/{}_logfile".format(self.out_dir, self.sort_name), "w") as f:
            f.write("\nargparse:{}\n".format(self.args))
            f.write("\nrunning directory:{}\n".format(self.out_dir))
            f.write("\nrunning command line:maf-swap {}/03_lastal/{}.maf |last-split |maf-swap |maf-sort > {}/04_sort/{}.maf\n".format(self.out_dir, self.sort_name, self.out_dir, self.sort_name))
        f.close()

        sp_sort_1_cmd = "maf-swap {}/03_lastal/{}.maf".format(self.out_dir, self.sort_name)
        sp_sort_2_cmd = "last-split"
        sp_sort_3_cmd = "maf-swap"
        sp_sort_4_cmd = "maf-sort"

        self.sp_sort_1 = sp.Popen(shlex.split(sp_sort_1_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
        self.sp_sort_2 = sp.Popen(shlex.split(sp_sort_2_cmd), stdin=self.sp_sort_1.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
        self.sp_sort_3 = sp.Popen(shlex.split(sp_sort_3_cmd), stdin=self.sp_sort_2.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
        self.sp_sort_4 = sp.Popen(shlex.split(sp_sort_4_cmd), stdin=self.sp_sort_3.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
        sp_sort_4_output = self.sp_sort_4.communicate()[0]
        with open("{}/04_sort/{}.maf".format(self.out_dir, self.sort_name), "a", encoding='utf-8') as f:
            f.write(sp_sort_4_output.decode('UTF-8').strip()+"\n")
        f.close()
        self.sp_sort_1.stdout.close()
        self.sp_sort_2.stdout.close()
        self.sp_sort_3.stdout.close()
        self.sp_sort_4.stdout.close()

        # check sort_4 return code
        if self.sp_sort_4.returncode == 0:
            with open("{}/04_sort/{}_logfile".format(self.out_dir, self.sort_name), "a", encoding='utf-8') as f:
                f.write("\nmaf-sort successfully finished\n")
            f.close()
            return (0)
        else:
            with open("{}/04_sort/error_{}".format(self.out_dir, self.sort_name), "a", encoding='utf-8') as f:
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
            return(1)


    def fasta_download_run(self):
        with Pool(self.threads) as p:
            p.map(self.fasta_download_func, self.accession_list)


    def fasta_swap_run(self):
        with Pool(self.parallel) as p:
            p.map(self.fasta_swap_func, zip(self.fna_gz_list, self.fasta_swap_name_list))


    def lastdb_run(self):
        self.lastdb_func()


    def last_train_run(self):
        with Pool(self.parallel) as p:
            p.map(self.last_train_func, self.last_train_name_list)


    def lastal_run(self):
        with Pool(self.parallel) as p:
            p.map(self.lastal_func, self.lastal_name_list)


    def sort_run(self):
        with Pool(self.parallel) as p:
            p.map(self.sort_func, self.sort_name_list)


if __name__ == '__main__':
    wga_in_one_step()
