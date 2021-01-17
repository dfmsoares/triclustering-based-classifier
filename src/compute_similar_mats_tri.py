import identify_triclusters as idt
from tricluster import Tricluster
import process_data as pdata
import pandas as pd
import numpy as np
import sys
from pathlib import Path
import os

if __name__ == "__main__":
    
    datafile = sys.argv[1]
    file_trics = sys.argv[2]
    path_matr = sys.argv[3]
    target_var = sys.argv[4]
    n = int(sys.argv[5])
    cat_feats = list(map(str, sys.argv[6].strip('[]').split(',')))
    cont_feats = list(map(str, sys.argv[7].strip('[]').split(',')))


    data_in = pdata.load_data(datafile)

    filename = datafile.split("/")[-1]

    sps = pdata.computeNSnapshots(data_in, 3, target_var)

    y = list()
    values = list()
    cols = list()
    for p in sps.keys():
        tp = data_in[p]
        for snaps in sps[p]:
            l = list()
            l.append(p)
            for e in snaps:
                i = e[0]
                l += list(tp[i].values())[:-1]
            values.append(l)
            y.append(e[1])

    cols.append("REF")
    for ti in range(0,n):
        for f in list(data_in[1][ti].keys())[:-1]:
            cols.append(str(ti)+f)
        
    mats = pd.DataFrame(data=values,
                    columns=cols)

    n_feats = len(list(data_in[1][ti].keys())[:-1])
    ps_tr = list()
    for e in mats.values:
        p_tric = Tricluster(n, n_feats, 1)
        i = 0
        f = 0
        for v in e[1:]:
            p_tric.addValue("T-"+str(i), "S-"+str(f), "G-" + str(int(e[0])), v)
            if f == n_feats-1:
                f = 0
                i +=1
            else:
                f+=1
        ps_tr.append(p_tric)


    if file_trics.endswith(".txt"):
        triclusters = idt.explore_triclusters(file_trics)
        matrix, cols = idt.compute_distance_matrix(ps_tr, triclusters, cat_feats, cont_feats, corr=False)
        idt.write_matrix(matrix, path_matr + "/" + filename[:-4] + "_DistanceMatrix.csv", y, cols, target_var)

        matrix, cols = idt.compute_distance_matrix(ps_tr, triclusters, cat_feats, cont_feats, corr=True)
        idt.write_matrix(matrix, path_matr + "/" + filename[:-4] + "_CorrelationMatrix.csv", y, cols, target_var)
