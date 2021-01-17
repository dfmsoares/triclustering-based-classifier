"""
Diogo Soares

Preprocess ALS dataset for triCluster Algorithm
"""

import pandas as pd
from pprint import pprint
import math
import datetime as dt

def load_data(file_name):
    """
    Load data from csv file into a Pandas DataFrame
    
    """
    data = pd.read_csv(file_name, sep=',')
    data_dict = data.to_dict('index')
    final_dict = dict()

    id_paciente_glob = 0
    time_counter = 0
    for k in data_dict.keys():
        ref = data_dict[k]['REF'] 
        if ref != id_paciente_glob:
            id_paciente_glob = ref
            time_counter = 0
        
        del data_dict[k]['REF']
        if ref not in final_dict:
            final_dict [ref] = {time_counter: data_dict[k]}
        else:
            final_dict[ref][time_counter] = data_dict[k]
        
        time_counter += 1

    return final_dict


def write_tab_file(data_dict, file_name, feat_list, tps):
    with open(file_name, 'wt') as out:
        out.write("Total Times:\t" + str(len(data_dict[list(data_dict.keys())[0]].keys())) + '\n')
        out.write("Total Samples:\t" + str(len(feat_list)) + '\n')
        out.write("Total Genes:\t" + str(len(data_dict.keys())) + '\n')

        for i in range(0, tps):
                out.write("Time\t" + str(i)+ '\n')
                out.write("ID\tNAME\t")
                for l in range(0, len(feat_list)):
                    if l != len(feat_list)-1:
                        out.write("S-" + str(l) + '\t')
                    else: 
                        out.write("S-" + str(l) + '\n')

                p = 0
                for (patient, t) in data_dict.items():
                        t_string = "{}\t"*(len(feat_list)+2)
                        tupl = (p, "G-" + str(patient)) + tuple(map(lambda x: t[i][x], feat_list))
                        line = t_string.format(*tupl)
                        out.write(line + '\n')
                        p+=1
        
    out.close()


def computeNSnapshots(data, n=3, label="Class"):
    final = dict()
    for (p, t) in data.items():
        if len(t.keys()) >= n:
            fd = dict()
            for (key, val) in t.items():
                fd[key] = val
                final[p] = fd
    
    snaps = dict()
    for (p, ts) in final.items():
        for t in ts.keys():
            if t < len(ts.keys())-(n-1) and all(map(lambda c: c != "Y", [final[p][t+y][label] for y in range(0,n-1)])):
                if p not in snaps:
                    snaps[p] = list()
                snaps[p].append([(t+j, final[p][t+j][label]) for j in range(0,n)])
    return snaps


def getClass(data, target_name):
    f = dict()
    for (p, t) in data.items():
        for (_, vals) in t.items():
            if p not in f:
                f[p] = [vals[target_name]]
            else:
                f[p].append(vals[target_name])

    return f
