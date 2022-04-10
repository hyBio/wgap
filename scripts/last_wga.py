#!/home/huyan/miniconda3/bin/python3
# _*_ coding: utf-8 _*_
# @Time : 2022/4/9 22:08
# @Author : 胡琰
# @Version：V 0.1
# @File : last_wga.py
# @Site :


import sys
import shlex
import subprocess as sp

class last_wga(object):

    def __init__(self):
        sys.stdout.write("<last_wga usage>: python last_wga.py <optional parameter> <subprocess_argparse>\n"+
                         "optional parameters:lastdb\tlast_train\tlastal\tlast_split\tmaf_swap\tmaf_sort\trename_species\tmultiz\n")
        if len(sys.argv) <= 1:
            sys.stderr.write("\nerror:No optional parameter given to last_wga"+
                             "\nExiting\n")
            sys.exit(0)
        else:
            self.optional_parameter = sys.argv[1]
            # check optional parameter
            if self.optional_parameter not in ["lastdb", "last_train", "lastal", "last_split", "maf_swap", "maf_sort", "rename_species", "multiz"]:
                sys.stderr.write("\nerror:Incorrect optional parameter {} given to last_wga\n".format(sys.argv[1]))
            else:
                # check subprocess argparse
                if len(sys.argv) < 3:
                    sys.stderr.write("\nerror:No sub process given to {}".format(self.optional_parameter)+
                                     "\nExiting"+
                                     "\nsubprocess help message:\n")
                    self.sp = sp.Popen(shlex.split("{} --help".format(self.optional_parameter.replace("_","-"))), stdout=sp.PIPE, stderr=sp.PIPE)
                    for stdout_line in iter(self.sp.stdout.readline, b''):
                        sys.stdout.write(stdout_line.decode('utf-8'))
                    self.sp.communicate()
                    sys.exit(0)
                elif len(sys.argv) >= 3:
                    # get subprocess_argparse
                    self.subprocess_argparse = " ".join(sys.argv[2:])
                    # run subprocess
                    self.subprocess = eval(self.optional_parameter+"()")
                    self.subprocess.running(self.subprocess_argparse)



class lastdb(object):

    def __init__(self):
        pass

    def running(self,subprocess_argparse=None):
        self.subprocess_argparse = subprocess_argparse
        sys.stdout.write("\nrunning command line:lastdb "+self.subprocess_argparse+"\n")

        self.sp = sp.Popen(shlex.split("lastdb {}".format(self.subprocess_argparse)), stdout=sp.PIPE, stderr=sp.PIPE)
        # check lastdb return code
        if self.sp.returncode == 0:
            sys.stdout.write("\nlastdb successfully finished\n")
        else:
            sys.stderr.write("\nerror:lastdb failed\n"+"lastdb error message:\n")
            for stderr_line in iter(self.sp.stderr.readline, b''):
                sys.stderr.write(stderr_line.decode('utf-8'))
            self.sp.communicate()

            self.sp_help = sp.Popen(shlex.split("lastdb --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            for stdout_line in iter(self.sp_help.stdout.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            self.sp_help.communicate()
            sys.exit(0)



class last_train(object):

    def __init__(self):
        pass

    def running(self,subprocess_argparse=None):
        self.subprocess_argparse = subprocess_argparse
        sys.stdout.write("\nrunning command line:last-train "+self.subprocess_argparse+"\n")
        self.sp = sp.Popen(shlex.split("last-train {}".format(self.subprocess_argparse)), stdout=sp.PIPE, stderr=sp.PIPE)
        # check last_train return code
        if self.sp.returncode == 0:
            sys.stdout.write("\nlast-train successfully finished\n")
        else:
            sys.stderr.write("\nerror:last-train failed\n"+"last-train error message:\n")
            for stderr_line in iter(self.sp.stderr.readline, b''):
                sys.stderr.write(stderr_line.decode('utf-8'))
            self.sp.communicate()

            self.sp_help = sp.Popen(shlex.split("last-train --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            for stdout_line in iter(self.sp_help.stdout.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            self.sp_help.communicate()
            sys.exit(0)



class lastal(object):

    def __init__(self):
        pass

    def running(self,subprocess_argparse=None):
        self.subprocess_argparse = subprocess_argparse
        sys.stdout.write("\nrunning command line:lastal "+self.subprocess_argparse+"\n")
        self.sp = sp.Popen(shlex.split("lastal {}".format(self.subprocess_argparse)), stdout=sp.PIPE, stderr=sp.PIPE)
        # check lastal return code
        if self.sp.returncode == 0:
            sys.stdout.write("\nlastal successfully finished\n")
        else:
            sys.stderr.write("\nerror:lastal failed\n"+"lastal error message:\n")
            for stderr_line in iter(self.sp.stderr.readline, b''):
                sys.stderr.write(stderr_line.decode('utf-8'))
            self.sp.communicate()

            self.sp_help = sp.Popen(shlex.split("lastal --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            for stdout_line in iter(self.sp_help.stdout.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            self.sp_help.communicate()
            sys.exit(0)



class last_split(object):

    def __init__(self):
        pass

    def running(self,subprocess_argparse=None):

        self.subprocess_argparse = subprocess_argparse
        sys.stdout.write("\nrunning command line:last-split "+self.subprocess_argparse+"\n")
        self.sp = sp.Popen(shlex.split("last-split {}".format(self.subprocess_argparse)), stdout=sp.PIPE, stderr=sp.PIPE)
        # check last_split return code
        if self.sp.returncode == 0:
            sys.stdout.write("\nlast-split successfully finished\n")
        else:
            sys.stderr.write("\nerror:last-split failed\n"+"last-split error message:\n")
            for stderr_line in iter(self.sp.stderr.readline, b''):
                sys.stderr.write(stderr_line.decode('utf-8'))
            self.sp.communicate()

            self.sp_help = sp.Popen(shlex.split("last-split --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            for stdout_line in iter(self.sp_help.stdout.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            self.sp_help.communicate()
            sys.exit(0)



class maf_swap(object):

    def __init__(self):
        pass

    def running(self,subprocess_argparse=None):

        self.subprocess_argparse = subprocess_argparse
        sys.stdout.write("\nrunning command line:maf-swap "+self.subprocess_argparse+"\n")
        self.sp = sp.Popen(shlex.split("maf-swap {}".format(self.subprocess_argparse)), stdout=sp.PIPE, stderr=sp.PIPE)
        # check maf_swap return code
        if self.sp.returncode == 0:
            sys.stdout.write("\nmaf-swap successfully finished\n")
        else:
            sys.stderr.write("\nerror:maf-swap failed\n"+"maf-swap error message:\n")
            for stderr_line in iter(self.sp.stderr.readline, b''):
                sys.stderr.write(stderr_line.decode('utf-8'))
            self.sp.communicate()

            self.sp_help = sp.Popen(shlex.split("maf-swap --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            for stdout_line in iter(self.sp_help.stdout.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            self.sp_help.communicate()
            sys.exit(0)



class maf_sort(object):

    def __init__(self):
        pass

    def running(self,subprocess_argparse=None):

        self.subprocess_argparse = subprocess_argparse
        sys.stdout.write("\nrunning command line:maf-sort "+self.subprocess_argparse+"\n")
        self.sp = sp.Popen(shlex.split("maf-sort {}".format(self.subprocess_argparse)), stdout=sp.PIPE, stderr=sp.PIPE)
        # check maf_sort return code
        if self.sp.returncode == 0:
            sys.stdout.write("\nmaf-sort successfully finished\n")
        else:
            sys.stderr.write("\nerror:maf-sort failed\n"+"maf-sort error message:\n")
            for stderr_line in iter(self.sp.stderr.readline, b''):
                sys.stderr.write(stderr_line.decode('utf-8'))
            self.sp.communicate()

            self.sp_help = sp.Popen(shlex.split("maf-sort --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            for stdout_line in iter(self.sp_help.stdout.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            self.sp_help.communicate()
            sys.exit(0)



class rename_species(object):

    def __init__(self):
        pass

    def running(self,subprocess_argparse=None):

        self.subprocess_argparse = subprocess_argparse
        sys.stdout.write("\nrunning command line:rename-species "+self.subprocess_argparse+"\n")
        self.sp = sp.Popen(shlex.split("rename-species {}".format(self.subprocess_argparse)), stdout=sp.PIPE, stderr=sp.PIPE)
        # check rename_species return code
        if self.sp.returncode == 0:
            sys.stdout.write("\nrename-species successfully finished\n")
        else:
            sys.stderr.write("\nerror:rename-species failed\n"+"rename-species error message:\n")
            for stderr_line in iter(self.sp.stderr.readline, b''):
                sys.stderr.write(stderr_line.decode('utf-8'))
            self.sp.communicate()

            self.sp_help = sp.Popen(shlex.split("rename-species --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            for stdout_line in iter(self.sp_help.stdout.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            self.sp_help.communicate()
            sys.exit(0)



class multiz(object):

    def __init__(self):
        pass

    def running(self,subprocess_argparse=None):

        self.subprocess_argparse = subprocess_argparse
        sys.stdout.write("\nrunning command line:multiz "+self.subprocess_argparse+"\n")
        self.sp = sp.Popen(shlex.split("multiz {}".format(self.subprocess_argparse)), stdout=sp.PIPE, stderr=sp.PIPE)
        # check multiz return code
        if self.sp.returncode == 0:
            sys.stdout.write("\nmultiz successfully finished\n")
        else:
            sys.stderr.write("\nerror:multiz failed\n"+"multiz error message:\n")
            for stderr_line in iter(self.sp.stderr.readline, b''):
                sys.stderr.write(stderr_line.decode('utf-8'))
            self.sp.communicate()

            self.sp_help = sp.Popen(shlex.split("multiz --help"), stdout=sp.PIPE, stderr=sp.PIPE)
            for stdout_line in iter(self.sp_help.stdout.readline, b''):
                sys.stdout.write(stdout_line.decode('utf-8'))
            self.sp_help.communicate()
            sys.exit(0)



if __name__ == '__main__':
    main = last_wga()



