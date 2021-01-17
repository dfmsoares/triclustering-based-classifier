import pandas as pd
import sys
import process_data as pdata

def write_tab_file(df, file_name, features, tps):
    with open(file_name, 'wt') as out:
        out.write("Total Times:\t" + str(tps) + '\n')
        out.write("Total Samples:\t" + str(len(features)) + '\n')
        out.write("Total Genes:\t" + str(len(df)) + '\n')

        for i in range(0,tps):
            fs = list(map(lambda x: str(i)+x, features))
            sub_feat = df[["REF"]+fs]
         
            out.write("Time\t" + str(i)+ '\n')
            out.write("ID\tNAME\t")
            for l in range(len(features)):
                if l != len(features)-1:
                    out.write("S-" + str(l) + '\t')
                else: 
                    out.write("S-" + str(l) + '\n')

            for index, row in sub_feat.iterrows():

                t_string = "{}\t"*(len(fs)+2)
                tupl = (index, "G-" + str(index)) + tuple(map(lambda x: row[x], fs))
                line = t_string.format(*tupl)
                out.write(line + '\n')

    out.close()

in_file = sys.argv[1]
out_file = sys.argv[2]
target_var = sys.argv[3]
n = int(sys.argv[4])

data_in = pdata.load_data(in_file)
sps = pdata.computeNSnapshots(data_in, n, target_var)

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


write_tab_file(mats, out_file, list(data_in[1][ti].keys())[:-1], n)