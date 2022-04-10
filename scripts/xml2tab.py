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
