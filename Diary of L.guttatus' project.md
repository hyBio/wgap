# Diary of *L.guttatus'* project

[TOC]



## 2022年3月21日15点06分

### 下载NCBI分类数据库

- 物种及分类信息文件：[https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz](https://links.jianshu.com/go?to=https%3A%2F%2Fftp.ncbi.nlm.nih.gov%2Fpub%2Ftaxonomy%2Ftaxdump.tar.gz)
- 解压文件 taxdump.tar.gz：

```shell
cd /data/huyan/nr_db
mkdir taxdump
cd taxdump/
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf ./taxdump.tar.gz
ll

total 453804
-rw-r--r-- 1 9019  583  19149208 Mar 21 15:28 citations.dmp
-rw-r--r-- 1 9019  583   4283609 Mar 21 15:26 delnodes.dmp
-rw-r--r-- 1 9019  583       452 Mar 21 15:20 division.dmp
-rw-r--r-- 1 9019  583     16444 Mar 21 15:28 gc.prt
-rw-r--r-- 1 9019  583      4921 Mar 21 15:20 gencode.dmp
-rw-r--r-- 1 9019  583   1239555 Mar 21 15:26 merged.dmp
-rw-r--r-- 1 9019  583 216720982 Mar 21 15:28 names.dmp
-rw-r--r-- 1 9019  583 165641529 Mar 21 15:27 nodes.dmp
-rw-rw---- 1 4544  583      2666 Sep 12  2019 readme.txt
-rw-r--r-- 1 root root  57612439 Mar 21 15:28 taxdump.tar.gz
```

- taxdump 目录中有两个重要文件：

  names.dmp：记录物种名及其分类编号

  nodes.dmp：记录分类编号的分类节点信息

### 功能注释

对于NR数据库，一般都是blast直接比对，获取基因或者蛋白质的功能注释信息，更近一步，如果有Taxonomy 数据库的数据，我们还可以做物种注释。

------

blast输出格式：（首选5）

-outfmt <String>
alignment view options:
0 = pairwise,显示具体匹配信息
1 = query-anchored showing identities,查询-比上区域，显示一致性
2 = query-anchored no identities,查询-比上区域，不显示一致性
3 = flat query-anchored, show identities,查询-比上区域的屏文形式，显示一致性
4 = flat query-anchored, no identities,查询-比上区域的屏文形式，不显示一致性
5 = XML Blast output,XML格式的输出
6 = tabular,TAB格式的输出
7 = tabular with comment lines,带注释行的TAB格式的输出
8 = Text ASN.1,文本方式的ASN格式输出
9 = Binary ASN.1,二进制方式的ASN格式输出  

------



```shell
screen -S lg

(base) [root@master l.guttatus]# pwd
/data/huyan/l.guttatus


#将fasta文件中多余的.删除
sed 's/\.$//g' /new_data/huyan/fish/lampris_guttatus/assembly/Chr_genome_final_gene.gff3.pep > ./Chr_genome_final_gene.gff3.pep

####
blastp -num_threads 10 -query ./Chr_genome_final_gene.gff3.pep -db /new_data/public_db/nr_blast_db/nr-2021-5-12/nr -evalue 1e-5 -num_alignments 5 -outfmt 5 -out /data/huyan/l.guttatus/l.guttatus.blastp.xml &
# 
[2] 38441
```

**速度太慢，弃用，改用diamond!**

```shell

diamond blastp -p 40 -d /new_data/public_db/nr_prot.dmnd -q ./Chr_genome_final_gene.gff3.pep -f 5 -o l.guttatus.blastp.xml &
[1] 120602
```



#### 格式转换（xml to txt）

```python
# _*_ coding: utf-8 _*_
# @Time : 2022/3/22 15:41 
# @Author : 胡琰 
# @Version：V 0.1
# @File : xml2tab.py
# @Site :


import xmltodict
import pandas as pd
import sys

xml_file = sys.argv[1]

def run():
    total_info=pd.DataFrame()
    with open(xml_file) as fd:
        doc = xmltodict.parse(fd.read())
        iterations = doc["BlastOutput"]["BlastOutput_iterations"]["Iteration"]

    for i in range(0,len(iterations)):
        iteration = iterations[i]

        if iteration["Iteration_hits"] is not None:
            if isinstance(iteration["Iteration_hits"]["Hit"],list):
                iteration_info = pd.concat([pd.DataFrame(iteration["Iteration_hits"]["Hit"][j]["Hit_hsps"]) for j in range(len(iteration["Iteration_hits"]["Hit"]))],axis=1).T
                for hit_item in ["Hit_num","Hit_id","Hit_def","Hit_accession","Hit_len"]:
                    iteration_info[hit_item] = [iteration["Iteration_hits"]["Hit"][j][hit_item] for j in range(len(iteration["Iteration_hits"]["Hit"]))]

            else:
                iteration_info = pd.DataFrame(iteration["Iteration_hits"]["Hit"]["Hit_hsps"]).T
                for hit_item in ["Hit_num","Hit_id","Hit_def","Hit_accession","Hit_len"]:
                    iteration_info[hit_item] = [iteration["Iteration_hits"]["Hit"][hit_item]]

            iteration_info["Iteration_iter_num"] = iteration["Iteration_iter-num"]
            iteration_info["Iteration_query_ID"] = iteration["Iteration_query-ID"]
            iteration_info["Iteration_query_def"] = iteration["Iteration_query-def"]
            iteration_info["Iteration_query_len"] = iteration["Iteration_query-len"]

            for key,values in iteration["Iteration_stat"]["Statistics"].items():
                iteration_info[key] = values

            iteration_info = iteration_info.drop(["Hsp_hseq","Hsp_qseq","Hsp_midline"],axis = 1)

            total_info = pd.concat([total_info, iteration_info], ignore_index=True)
            iteration_info.to_csv("tmp.txt", sep='\t', index=None, header=False, mode="a")

        else:
            pass

    total_info.to_csv("lg_xml.txt",sep='\t',index=None)


if __name__ == '__main__':
    run()
```

```shell
[huyan@master lg_gene]$ python3 xml2tab.py "/data/huyan/l.guttatus/l.guttatus.blastp.xml" &
[1] 109988
```



## 2022年3月28日10点46分

### 下载11个目的鱼的基因组信息

先准备好11个目的鱼的信息，如下：

| Order              | Organism Name              | Assembly Level | Assembly  Accession | Assembly Name                      | Assembly Stats  Total Sequence Length | Assembly  Submission Date |
| ------------------ | -------------------------- | -------------- | ------------------- | ---------------------------------- | ------------------------------------- | ------------------------- |
| Zeiformes          | Cyttopsis rosea            | Scaffold       | GCA_900302355.1     | ASM90030235v1                      | 546506150                             | 2018/3/19                 |
| Zeiformes          | Zeus faber                 | Scaffold       | GCA_900323335.1     | ASM90032333v1                      | 610433400                             | 2018/3/26                 |
| Trachichthyiformes | Anoplogaster cornuta       | Scaffold       | GCA_900683385.1     | Anoplogaster_cornuta_assembly      | 404844788                             | 2019/6/15                 |
| Trachichthyiformes | Diretmus argenteus         | Scaffold       | GCA_900660295.1     | Diretmus_argenteus_assembly        | 302363458                             | 2019/6/15                 |
| Trachichthyiformes | Hoplostethus  atlanticus   | Scaffold       | GCA_900660355.1     | Hoplostethus_atlanticus_assembly   | 520173038                             | 2019/6/18                 |
| Trachichthyiformes | Monocentris japonicus      | Scaffold       | GCA_900323365.1     | ASM90032336v1                      | 556023515                             | 2018/3/26                 |
| Stylephoriformes   | Stylephorus chordatus      | Scaffold       | GCA_900312615.1     | ASM90031261v1                      | 488488587                             | 2018/3/20                 |
| Scombriformes      | Euthynnus affinis          | Scaffold       | GCA_019973915.1     | EUTAFF_01                          | 758243246                             | 2021/9/9                  |
| Scombriformes      | Pampus argenteus           | Scaffold       | GCA_000697985.1     | PamArg1.0                          | 350448509                             | 2014/6/3                  |
| Scombriformes      | Scomber colias             | Scaffold       | GCA_021039105.1     | S.colias_p1.0                      | 814072382                             | 2021/12/3                 |
| Scombriformes      | Thunnus albacares          | Chromosome     | GCF_914725855.1     | fThuAlb1.1                         | 792101331                             | 2021/9/16                 |
| Scombriformes      | Thunnus albacares          | Contig         | GCA_914744365.1     | fThuAlb1.1 alternate  haplotype    | 786643039                             | 2021/9/15                 |
| Scombriformes      | Thunnus maccoyii           | Chromosome     | GCF_910596095.1     | fThuMac1.1                         | 782407452                             | 2021/7/7                  |
| Scombriformes      | Thunnus maccoyii           | Contig         | GCA_910596085.1     | fThuMac1.1 alternate  haplotype    | 782574410                             | 2021/7/7                  |
| Polymixiiformes    | Polymixia japonica         | Scaffold       | GCA_900302305.1     | ASM90030230v1                      | 554895936                             | 2018/3/19                 |
| Ophidiiformes      | Brotula barbata            | Scaffold       | GCA_900303265.1     | ASM90030326v1                      | 485061200                             | 2018/3/20                 |
| Ophidiiformes      | Carapus acus               | Scaffold       | GCA_900312935.1     | ASM90031293v1                      | 387834307                             | 2018/3/18                 |
| Ophidiiformes      | Lamprogrammus exutus       | Scaffold       | GCA_900312555.1     | ASM90031255v1                      | 492850272                             | 2018/3/20                 |
| Ophidiiformes      | Lucifuga dentata           | Scaffold       | GCA_014773175.1     | Ldentata1.0                        | 630651555                             | 2020/9/28                 |
| Lampriformes       | Lampris guttatus           | Scaffold       | GCA_900302545.1     | ASM90030254v1                      | 849277706                             | 2018/3/30                 |
| Lampriformes       | Lampris incognitus         | Chromosome     | GCA_022059245.1     | ASM2205924v1                       | 1367434603                            | 2022/2/10                 |
| Lampriformes       | Lampris megalopsis         | Contig         | GCA_022114975.1     | ASM2211497v1                       | 1090846196                            | 2022/2/9                  |
| Lampriformes       | Regalecus glesne           | Scaffold       | GCA_900302585.1     | ASM90030258v1                      | 656003707                             | 2018/3/20                 |
| Lamniformes        | Carcharodon  carcharias    | Chromosome     | GCF_017639515.1     | sCarCar2.pri                       | 4286294447                            | 2021/3/30                 |
| Holocentriformes   | Holocentrus rufus          | Scaffold       | GCA_900302615.1     | ASM90030261v1                      | 649757301                             | 2018/3/19                 |
| Holocentriformes   | Myripristis jacobus        | Scaffold       | GCA_900302555.1     | ASM90030255v1                      | 720396841                             | 2018/3/21                 |
| Holocentriformes   | Myripristis murdjan        | Chromosome     | GCF_902150065.1     | fMyrMur1.1                         | 835254674                             | 2019/7/3                  |
| Gadiformes         | Bathygadus  melanobranchus | Scaffold       | GCA_900302375.1     | ASM90030237v1                      | 431202967                             | 2018/3/19                 |
| Gadiformes         | Gadiculus argenteus        | Scaffold       | GCA_900302595.1     | ASM90030259v1                      | 396767394                             | 2018/3/18                 |
| Gadiformes         | Gadus chalcogrammus        | Scaffold       | GCA_900302575.1     | ASM90030257v1                      | 448868398                             | 2018/3/19                 |
| Gadiformes         | Gadus morhua               | Chromosome     | GCF_902167405.1     | gadMor3.0                          | 669949713                             | 2019/7/16                 |
| Gadiformes         | Macrourus berglax          | Scaffold       | GCA_900302365.1     | ASM90030236v1                      | 399875629                             | 2018/3/18                 |
| Gadiformes         | Merlangius merlangus       | Scaffold       | GCA_900323355.1     | ASM90032335v1                      | 423942190                             | 2018/3/26                 |
| Gadiformes         | Merluccius capensis        | Scaffold       | GCA_900312945.1     | ASM90031294v1                      | 414317329                             | 2018/3/18                 |
| Gadiformes         | Merluccius  merluccius     | Scaffold       | GCA_900312545.1     | ASM90031254v1                      | 401034705                             | 2018/3/18                 |
| Gadiformes         | Phycis blennoides          | Scaffold       | GCA_900302315.1     | ASM90030231v1                      | 416766999                             | 2018/3/18                 |
| Gadiformes         | Phycis phycis              | Scaffold       | GCA_900302335.1     | ASM90030233v1                      | 346335180                             | 2018/3/18                 |
| Cetomimiformes     | Rondeletia loricata        | Scaffold       | GCA_900302605.1     | ASM90030260v1                      | 568597941                             | 2018/3/19                 |
| Beryciformes       | Diretmoides  pauciradiatus | Scaffold       | GCA_900660315.1     | Diretmoides_pauciradiatus_assembly | 672603014                             | 2019/6/19                 |
| Beryciformes       | Gephyroberyx  darwinii     | Scaffold       | GCA_900660455.1     | Gephyroberyx_darwinii_assembly     | 535045983                             | 2019/6/18                 |
| Beryciformes       | Acanthochaenus  luetkenii  | Scaffold       | GCA_900312575.1     | ASM90031257v1                      | 545759480                             | 2018/3/19                 |
| Beryciformes       | Beryx splendens            | Scaffold       | GCA_900312565.1     | ASM90031256v1                      | 533267752                             | 2018/3/19                 |
| Beryciformes       | Neoniphon sammara          | Scaffold       | GCA_900302535.1     | ASM90030253v1                      | 659225728                             | 2018/3/21                 |
| Batrachoidiformes  | Chatrabus melanurus        | Scaffold       | GCA_900302635.1     | ASM90030263v1                      | 1128442880                            | 2018/3/30                 |
| Batrachoidiformes  | Opsanus beta               | Scaffold       | GCA_900660325.1     | Opsanus_beta_assembly              | 1028783780                            | 2019/6/19                 |
| Batrachoidiformes  | Thalassophryne  amazonica  | Chromosome     | GCF_902500255.1     | fThaAma1.1                         | 2446592878                            | 2019/9/30                 |
| Batrachoidiformes  | Thalassophryne  amazonica  | Scaffold       | GCA_902500245.1     | fThaAma1.1 alternate  haplotype    | 1767810094                            | 2019/9/30                 |

将Assembly Accession整理到acc.list表格中，使用rsync进行批量安装，acc.list如下：

```SHELL
[huyan@master lampris]$ more acc.list
GCA_900302355.1
GCA_900323335.1
GCA_900683385.1
GCA_900660295.1
GCA_900660355.1
GCA_900323365.1
GCA_900312615.1
...
```

命令行如下：

```shell
for sample in `cat "/data/huyan/lampris/acc.list"`
do
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/${sample:0:3}/${sample:4:3}/${sample:7:3}/${sample:10:3}/${sample}* /data/huyan/lampris/ &
done
```

### 全基因组树的构建

#### 预处理

##### 物种名称缩写替换

```shell
# 为了防止合并基因组两两比对后，结果的易看性，先对物种基因组fa/gff文件中的序列名称前加上物种名称缩写，并删除线粒体序列
## 更改序列名称示例，GCF_914725855.1_fThuAlb1.1_genomic.fna.gz为使用的基因组文件名，缩写一般使用的属名和种名的前3个字母。
### 更改fasta文件中每个染色体的编号，加上物种的前缀，删掉多余的信息

for s in `cat "/data/huyan/lampris/fna.list"`
do
gca=(${s//// })
gca=${gca[3]}
gca=(${gca//./ })
gca=${gca[0]}.1
echo $gca

name=`cat "/data/huyan/lampris/acc2name.list"|grep $gca`
name=(${name// / })
name=${name[2]}
echo $name
zcat ${s}|awk -v FS=' ' '{print $1}' |sed "s/>/>${name}_/g" > ${name}_genome.fa &
done
```

##### 线粒体基因组提取

```shell
## 提取线粒体基因组序列,此步骤可根据实际情况进行省略
grep -A 300 ">SpaAur_NC_024236" fSpaAur_genome.fa  > fSpaAur_mito.fa
```

##### fasta文件建立索引

```shell
for sample in fCytFab fZeuFab fAnoCor fDirArg fHopAtl fMonJap fStyCho fEutAff fPamArg fScoCol fThuAlb1 fThuAlb2 fThuMac1 fThuMac2 fPolJap fBroBar fCarAcu fLamExu fLucDen fLamGut fLamInc fLamMeg fRegGle fCarCar fHolRuf fMyrJac fMyrMur fBatMel fGadArg fGadCha fGadMor fMacBer fMerMerla fMerCap fMerMerlu fPhyBle fPhyPhy fRonLor fDirPau fGepDar fAcaLue fBerSpl fNeoSam fChaMel fOpsBet fThaAma1 fThaAma2
do
samtools faidx ${sample}_genome.fa &
done
```

#### 建库

##### 创建参考物种基因组index

```shell
# 创建参考物种基因组index
lastdb -P20 -uNEAR -cR11 fThuMac1_db "/data/huyan/lampris/align_result/fThuMac1_genome.fa" &
```

[1] 4753

##### 参数训练（可略）

```shell
#Try to find suitable score parameters for aligning the given sequences.

for sample in fCytFab fZeuFab fAnoCor fDirArg fHopAtl fMonJap fStyCho fEutAff fPamArg fScoCol fThuAlb1 fThuAlb2 fThuMac1 fThuMac2 fPolJap fBroBar fCarAcu fLamExu fLucDen fLamGut fLamInc fLamMeg fRegGle fCarCar fHolRuf fMyrJac fMyrMur fBatMel fGadArg fGadCha fGadMor fMacBer fMerMerla fMerCap fMerMerlu fPhyBle fPhyPhy fRonLor fDirPau fGepDar fAcaLue fBerSpl fNeoSam fChaMel fOpsBet fThaAma1 fThaAma2
do
last-train -P 2 --revsym --matsym --gapsym -E0.05 -C2 /data/huyan/lampris/wga_db/fThuMac1_db /data/huyan/lampris/align_result/${sample}_genome.fa > /data/huyan/lampris/last_train/${sample}.mat &
done
```

[1] 73898 [2] 73899 [3] 73900 [4] 73901 [5] 73902 [6] 73903 [7] 73904 [8] 73905 [9] 73906 [10] 73907 [11] 73908 [12] 73910 [13] 73911 [14] 73912 [15] 73913 [16] 73914 [17] 73915 [18] 73916 [19] 73917 [20] 73918 [21] 73919 [22] 73920 [23] 73921 [24] 73922 [25] 73923 [26] 73924 [27] 73925 [28] 73926 [29] 73927 [30] 73928 [31] 73929 [32] 73930 [33] 73931 [34] 73932 [35] 73933 [36] 73934 [37] 73935 [38] 73936 [39] 73937 [40] 73938 [41] 73939 [42] 73940 [43] 73941 [44] 73943 [45] 73944 [46] 73945 [47] 73946 

##### 比对

```shell
cd /data/huyan/lampris/lastal
for sample in fCytFab fZeuFab fAnoCor fDirArg fHopAtl fMonJap fStyCho fEutAff fPamArg fScoCol fThuAlb1 fThuAlb2 fThuMac1 fThuMac2 fPolJap fBroBar fCarAcu fLamExu fLucDen fLamGut fLamInc fLamMeg fRegGle fCarCar fHolRuf fMyrJac fMyrMur fBatMel fGadArg fGadCha fGadMor fMacBer fMerMerla fMerCap fMerMerlu fPhyBle fPhyPhy fRonLor fDirPau fGepDar fAcaLue fBerSpl fNeoSam fChaMel fOpsBet fThaAma1 fThaAma2
do
lastal -m100 -E0.05 -C2 -P 40 -p /data/huyan/lampris/last_train/${sample}.mat /data/huyan/lampris/wga_db/fThuMac1_db /data/huyan/lampris/align_result/${sample}_genome.fa | last-split -m1 >${sample}.maf
done &
```























