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
        self.parser.add_argument('-b', '--begin', default="fasta_swap", help='the begin subprocess of the pipeline\n'+
                                                                             'default: fasta_swap\n'+
                                                                             'optional parameters:fasta_download\tfasta_swap\tlastdb\tlast_train\tlastal\tmultiz')
        self.parser.add_argument('-t', '--threads', type=int, default=1, help='number of threads')
        self.parser.add_argument('-c', '--configure', help='configure file path')
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
        self.begin = args.begin
        self.threads = self.args.threads
        self.configure = self.args.configure
        self.fa_dir = self.args.fa_dir
        self.ref_fa = self.args.ref_fa
        self.exclude_fa = self.args.exclude_fa
        self.out_dir = self.args.out_dir
        self.lastdb = self.args.lastdb
        self.last_train = self.args.last_train
        self.lastal = self.args.lastal

        if self.begin == "fasta_download":
            self.fasta_download()
            self.fasta_swap()
            self.lastdb_func()
            self.last_train_func()
            self.lastal_func()
            self.multiz_func()
        elif self.begin == "fasta_swap":
            self.fasta_swap()
            self.lastdb_func()
            self.last_train_func()
            self.lastal_func()
            self.multiz_func()
        elif self.begin == "lastdb":
            self.lastdb_func()
            # self.last_train_func()
            # self.lastal_func()
            # self.multiz_func()
        elif self.begin == "last_train":
            self.last_train_func()
            # self.lastal_func()
            # self.multiz_func()
        elif self.begin == "lastal":
            self.lastal_func()
            # self.multiz_func()
        elif self.begin == "multiz":
            self.multiz_func()
        else:
            print("Incorrect parameter '{}' given to whole_genome_alignmen, please check the help".format(self.begin))


    def fasta_swap(self):
        sp.Popen(shlex.split("mkdir -p {}/00_assembly_fasta".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/00_assembly_fasta/".format(self.out_dir)))
        os.chdir(r"{}00_assembly_fasta/".format(self.out_dir))
        with open(self.configure, 'r') as f:
            lines = f.readlines()
        for line in lines:
            self.fna_gz = line.split('\t')[0].strip()
            self.name = line.split('\t')[1].strip()

            sys.stdout = open("{}/00_assembly_fasta/{}_logfile".format(self.out_dir,self.name), "w", encoding='utf-8')
            sys.stdout.write("\nargparse:{}\n".format(self.args))
            sys.stdout.write("\nrunning directory:{}\n".format(self.out_dir))
            sys.stdout.write(" ".join(['\nrunning command line:', 'zcat', '{}'.format(self.fna_gz),'|', 'awk', '\'{print $1}\'', '|', 'sed', '\'s/>/>{}_/g\''.format(self.name), '>', '{}/00_assembly_fasta/{}.fa\n'.format(self.out_dir,self.name)]))
            sys.stdout.flush()

            self.sp_ungzip = sp.Popen(["zcat", "{}".format(self.fna_gz)], stdout=sp.PIPE)
            self.sp_awk = sp.Popen(shlex.split("awk '{print $1}'"), stdin=self.sp_ungzip.stdout, stdout=sp.PIPE)
            self.sp_ungzip.stdout.close()
            self.sp_sed = sp.Popen(shlex.split("sed 's/>/>{}_/g'".format(self.name)), stdin=self.sp_awk.stdout, stdout=sp.PIPE)
            self.sp_awk.stdout.close()
            fasta_swap_1_output = self.sp_sed.communicate()[0]
            with open('{}/00_assembly_fasta/{}.fa'.format(self.out_dir,self.name), 'a', encoding='utf-8') as fa:
                fa.write(fasta_swap_1_output.decode('UTF-8').strip() + '\n')
            self.sp_sed.stdout.close()

            if self.sp_sed.returncode == 0:
                sys.stdout.write("\nfasta_swap_1 successfully finished\n")
                sys.stdout.write("\nrunning command line:samtools faidx -@ {} {}/00_assembly_fasta/{}.fa\n".format(self.threads, self.out_dir, self.name))
                sys.stdout.flush()
                self.sp_faidx = sp.Popen(shlex.split("samtools faidx {}/00_assembly_fasta/{}.fa".format(self.out_dir,self.name)), stdout=sp.PIPE)
                fasta_swap_2_output = self.sp_faidx.communicate()[0]
                with open('{}/00_assembly_fasta/{}.fa.fai'.format(self.out_dir,self.name), 'a', encoding='utf-8') as fai:
                    fai.write(fasta_swap_2_output.decode('UTF-8').strip() + '\n')
                self.sp_faidx.stdout.close()

                if self.sp_faidx.returncode == 0:
                    sys.stdout.write("\nfasta_swap_2 successfully finished\n")
                    sys.stdout.flush()
                    sys.stdout.close()
                    return(0)
                else:
                    sys.stdout.write("\nfasta_swap_2 failed\n")
                    sys.stdout.flush()
                    sys.stdout.close()
                    return(1)
            else:
                sys.stdout.write("\nfasta_swap_1 failed\n")
                sys.stdout.flush()
                sys.stdout.close()
                return(1)


    def lastdb_func(self):
        sp.Popen(shlex.split("mkdir -p {}/01_lastdb".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/01_lastdb/".format(self.out_dir)))

        sys.stdout = open("{}/01_lastdb/logfile".format(self.out_dir), "w", encoding='utf-8')
        sys.stdout.write("\nargparse:{}\n".format(self.args))
        sys.stdout.write("\nrunning directory:{}\n".format(self.out_dir))
        sys.stdout.write("\nrunning command line:lastdb -P {} {} {}_db {}/00_assembly_fasta/{}.fa\n".format(self.threads, self.lastdb, self.ref_fa, self.out_dir, self.ref_fa))
        sys.stdout.flush()

        os.chdir(r"{}/01_lastdb/".format(self.out_dir))
        sp_lastdb_cmd = "lastdb -P {} ".format(self.threads) + self.lastdb + " " + self.ref_fa + '_db ' + self.out_dir + '/00_assembly_fasta/' + self.ref_fa + '.fa'
        self.sp_lastdb = sp.Popen(shlex.split(sp_lastdb_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
        # wiat for the process to finish
        self.sp_lastdb.communicate()
        # check lastdb return code
        if self.sp_lastdb.returncode == 0:
            sys.stdout = open("{}/01_lastdb/logfile".format(self.out_dir), "a", encoding='utf-8')
            sys.stdout.write("\nlastdb successfully finished\n")
            sys.stdout.flush()
            return(0)
        else:
            sys.stdout = open("{}/01_lastdb/logfile".format(self.out_dir), "a", encoding='utf-8')
            sys.stdout.write("\nerror:lastdb failed\n"+"\nlastdb error message:\n")
            sys.stdout.flush()

            sys.stderr = open("{}/01_lastdb/logfile".format(self.out_dir), "a", encoding='utf-8')
            sys.stderr.write(self.sp_lastdb.communicate()[1].decode('UTF-8'))
            sys.stderr.flush()

            self.sp_lastdb_help = sp.Popen(shlex.split("lastdb --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            with open("{}/01_lastdb/logfile".format(self.out_dir), "a") as f:
                f.write(self.sp_lastdb_help.stdout.read().decode('UTF-8'))
                f.write(self.sp_lastdb_help.stderr.read().decode('UTF-8'))
            self.sp_lastdb_help.stdout.close()
            self.sp_lastdb_help.stderr.close()
            return(1)


    def last_train_func(self):
        sp.Popen(shlex.split("mkdir -p {}/02_last_train".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/02_last_train/".format(self.out_dir)))
        os.chdir(r"{}/02_last_train/".format(self.out_dir))

        with open(self.configure, 'r') as f:
            lines = f.readlines()
        for line in lines:
            self.name = line.split('\t')[1].strip()
            if self.name != self.ref_fa:
                sys.stdout = open("{}/02_last_train/{}_logfile".format(self.out_dir,self.name), "w", encoding='utf-8')
                sys.stdout.write("\nargparse:{}\n".format(self.args))
                sys.stdout.write("\nrunning directory:{}\n".format(self.out_dir))
                sys.stdout.write("\nrunning command line:last-train -P {} {} {}/01_lastdb/{}_db {}/00_assembly_fasta/{}.fa > {}/02_last_train/{}.mat\n".format(self.threads, self.last_train, self.out_dir, self.ref_fa, self.out_dir, self.name, self.out_dir, self.name))
                sys.stdout.flush()

                sp_last_train_cmd = "last-train -P {} ".format(self.threads) + self.last_train + " " + self.out_dir + "/01_lastdb/" + self.ref_fa + "_db " + self.out_dir + "/00_assembly_fasta/" + self.name + ".fa"
                self.sp_last_train = sp.Popen(shlex.split(sp_last_train_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
                # wiat for the process to finish
                last_train_output = self.sp_last_train.communicate()[0]
                with open("{}/02_last_train/{}.mat".format(self.out_dir, self.name), "a") as f:
                    f.write(last_train_output.decode('UTF-8').strip()+"\n")
                self.sp_last_train.stdout.close()

                # check lastdb return code
                if self.sp_last_train.returncode == 0:
                    sys.stdout = open("{}/02_last_train/{}_logfile".format(self.out_dir,self.name), "a", encoding='utf-8')
                    sys.stdout.write("\nlast-train successfully finished\n")
                    sys.stdout.flush()
                else:
                    sys.stderr = open("{}/02_last_train/error_{}".format(self.out_dir,self.name), "a", encoding='utf-8')
                    sys.stderr.write("\nerror:last-train failed\n"+"last-train error message:\n")
                    sys.stderr.write(self.sp_last_train.communicate()[1].decode('UTF-8'))
                    sys.stderr.flush()


    def lastal_func(self):
        sp.Popen(shlex.split("mkdir -p {}/03_lastal".format(self.out_dir)))
        sp.Popen(shlex.split("cd {}/03_lastal/".format(self.out_dir)))
        os.chdir(r"{}/03_lastal/".format(self.out_dir))

        with open(self.configure, 'r') as f:
            lines = f.readlines()
        for line in lines:
            self.name = line.split('\t')[1].strip()
            if self.name != self.ref_fa:
                with open("{}/03_lastal/{}_logfile".format(self.out_dir, self.name), "w") as f:
                    f.write("\nargparse:{}\n".format(self.args))
                    f.write("\nrunning directory:{}\n".format(self.out_dir))
                    f.write("\nrunning command line:lastal -P {} {} {}/02_last_train/{}.mat {}/01_lastdb/{}_db {}/00_assembly_fasta/{}.fa > {}/03_lastal/{}.maf\n".format(self.threads, self.lastal, self.out_dir, self.name, self.out_dir, self.name, self.out_dir, self.name, self.out_dir, self.name))
                f.close()

                sp_lastal_1_cmd = "lastal -P {} ".format(self.threads) + self.lastal + " " + self.out_dir + "/02_last_train/" + self.name + ".mat " + self.out_dir + "/01_lastdb/" + self.ref_fa + "_db " + self.out_dir + "/00_assembly_fasta/" + self.name + ".fa "
                sp_lastal_2_cmd ="last-split " + "-f " + "MAF+"

                self.sp_lastal_1 = sp.Popen(shlex.split(sp_lastal_1_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
                self.sp_lastal_1.communicate()
                # check lastal_1 return code
                if self.sp_lastal_1.returncode == 0:
                    with open("{}/03_lastal/{}_logfile".format(self.out_dir, self.name), "a", encoding='utf-8') as f:
                        f.write("\nlastal successfully finished\n")
                    f.close()
                else:
                    with open("{}/03_lastal/{}_logfile".format(self.out_dir, self.name), "a", encoding='utf-8') as f:
                        f.write("\nerror:lastal failed\n" + "lastal_1 error message:\n")
                        f.write(self.sp_lastal_1.stderr.read().decode('UTF-8'))
                        self.sp_lastal_1.stderr.close()
                    f.close()
                    return(1)

                self.sp_lastal_2 = sp.Popen(shlex.split(sp_lastal_2_cmd), stdin=self.sp_lastal_1.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
                with open("{}/03_lastal/{}.maf".format(self.out_dir, self.name), "a") as f:
                    for line in self.sp_lastal_2.stdout:
                        f.write(line.decode('UTF-8'))
                f.close()
                self.sp_lastal_2.stdout.close()
                self.sp_lastal_1.stdout.close()

                # check lastal_2 return code
                if self.sp_lastal_2.returncode == 0:
                    with open("{}/03_lastal/{}_logfile".format(self.out_dir, self.name), "a", encoding='utf-8') as f:
                        f.write("\nlast-split successfully finished\n")
                    f.close()
                    return(0)
                else:
                    with open("{}/03_lastal/{}_logfile".format(self.out_dir, self.name), "a", encoding='utf-8') as f:
                        f.write("\nerror:last-split failed\n" + "last-split error message:\n")
                        f.write(self.sp_lastal_1.stderr.read().decode('UTF-8'))
                        self.sp_lastal_1.stderr.close()
                    f.close()
                    return(1)




if __name__ == '__main__':
    wga_in_one_step()
