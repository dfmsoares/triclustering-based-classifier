from tricluster import Tricluster
import itertools
import scipy.spatial.distance as distance
import scipy.stats as ss
import statistics as stat
import pandas as pd
import numpy as np
import math
import copy


def getOverlaping(output_file):
    fp = open(output_file, 'r')
    content = fp.readlines()

    content = list(map(lambda x: x.strip(), content))
    ok = False
    for c in content:
        if c.startswith("Overlap") and not ok:
            ok = True
        elif c.startswith("Overlap") and ok:
            over = c.split(":")[1]
            over = over.strip()
            over = float(over[:-1])
            return over


def explore_triclusters(output_file):
    fp = open(output_file, 'r')
    content = fp.readlines()

    content = list(map(lambda x: x.strip(), content))

    flag = False
    triclusters = list()
    for c in content:
        if c.startswith("====================================================================================================================="):
            flag = True
        elif c.startswith("|T|x|S|x|G|") and flag:
            l = c.split(":")[1].split('x')
            tri = Tricluster(int(l[0]), int(l[1]), int(l[2]))
        elif c.startswith("Time") and flag:
            time = int(c.split(":")[1].strip())
            tri.addTime("T-" + str(time))
        elif c.startswith("S-") and flag:
            smp = list(map(lambda x: int(x.strip()), c.split("S-")[1:]))
            for s in smp:
                tri.addSample("S-" + str(s))
        elif c.startswith("G") and flag:
            pc = c.split('\t')[0]
            pci = int(pc.split("-")[1])
            tri.addPatient(pc)
            q = c.split("\t")[1].strip().replace(" ", "$")
            vals = q.split("$")
            vals = list(filter(lambda z: len(z) > 0, vals))
            idx = 0
            for v in vals:
                s = tri.getSamples()[idx]
                tri.addValue("T-" + str(time), s, pc, float(v.strip()))
                idx += 1
        elif c.startswith("Cluster") and flag:
            triclusters.append(tri)
    return triclusters


def compute_representative_patterns(tricluster, categorical_feats, continuos_feats):
    """
    Format: [[], []]
    """
    def aux_mean(z):
        z = list(filter(lambda x: x != 0, z))
        if len(z) == 0:
            return 0
        else:
            return stat.mean(z)

    def aux_mode(z):
        if len(z) == 0:
            return 0
        else:
            return ss.mode(z, keepdims=True)[0][0]

    result = []
    for tp in tricluster.getTimes():
        bic_vals = list()
        bicluster = [tricluster.getFeatValues(
            tp, f) for f in tricluster.getSamples()]

        # print(bicluster)
        means = list(map(lambda z: aux_mean(z), bicluster))
        modes = list(map(lambda z: aux_mode(z), bicluster))

        # for l in bicluster:
        #     means.append(stat.mean(l))
        #     modes.append(ss.mode(l)[0][0])
        i = 0
        for ft in tricluster.getSamples():
            if ft in categorical_feats:
                bic_vals.append((ft, modes[i]))
            else:
                bic_vals.append((ft, means[i]))
            i += 1
        result.append(bic_vals)

    # return np.concatenate((mode, mean), axis=1)

    return (result)


def compute_matrix(patient_ids, triclusters):
    final_matrix = list()
    for p in patient_ids:
        line = [1 if p in t.getGPatients() else 0 for t in triclusters]
        final_matrix.append(line)

    cols = ["Tric_"+str(it) for it in range(len(triclusters))]

    return final_matrix, cols


def compute_MSR_matrix(patients_tri, triclusters):
    """
        patients_tri is a list with Tricluster elements, representing a Patient Slice of 
            original data
        triclusters is a list with Tricluster elements, representing a Tricluster
    """

    final_matrix = list()
    t_msr = [t.compute_MSR() for t in triclusters]
    p_msr = [p.compute_MSR() for p in patients_tri]
    for p in p_msr:
        line = [abs(p - t) for t in t_msr]
        final_matrix.append(line)

    return final_matrix


def compute_Correlation_matrix(patients_tri, triclusters):
    final_matrix = list()
    print(len(patients_tri))
    for bic_p in patients_tri:
        line = list()
        for tric in triclusters:
            corr_list = list()
            for tric_p in tric.getPatients():
                tric_tps = tric.getTimes()
                tric_fs = tric.getSamples()
                for fo in tric_fs:
                    p_slice = bic_p.getSlice(c=fo)[:len(tric_tps)]
                    # print("fop", p_slice)
                    t_slice = tric.getSlice(g=tric_p, c=fo)
                    # print("fot", t_slice)

                    corr_list.append(1 - distance.canberra(p_slice, t_slice))

                    # print(corr_list)
                for to in tric_tps:
                    p_slice = bic_p.getSlice(t=to)[:len(tric_fs)]
                    # print("top", p_slice)
                    t_slice = tric.getSlice(g=tric_p, t=to)
                    # print("tot", t_slice)
                    #
                    corr_list.append(1 - distance.canberra(p_slice, t_slice))
                    # print(corr_list)
            # rint(corr_list)
            line.append(stat.mean(corr_list))
        # print(len(line),line)
        final_matrix.append(line)

    return final_matrix


def compute_distance_matrix(patients_tri, triclusters, categorical_feats=[], continuos_feats=[], corr=False):
    dist_set = dict()
    print(len(patients_tri))

    p_pats_reps = list()
    for bic_p in patients_tri:
        # Compute Patient Patterns
        p_pats = [bic_p.getPatientsVals(t=tp) for tp in bic_p.getTimes()]

        p_pats_reps.append(p_pats)

    tric_names = list()
    t_pats_reps = list()
    i = 0
    for tric in triclusters:
        t_pats = compute_representative_patterns(
            tric, categorical_feats, continuos_feats)
        if t_pats not in t_pats_reps:
            t_pats_reps.append(t_pats)
            tric_names.append("Tric_" + str(i))
        i += 1
    pat = 0

    for p_pats in p_pats_reps:

        nTr = 0
        for t in t_pats_reps:

            ss = list(map(lambda x: int(x[0].split('-')[-1]), t[0]))

            for i in range(len(t)):
                a = np.array(tuple(list(map(lambda y: y[1], t[i]))))
                b = np.array(tuple(p_pats[i]))[[ss]]

                if not corr:
                    dist_set[(pat, tric_names[nTr], i)] = np.linalg.norm(a-b)
                else:
                    if len(a) == 2:
                        dist_set[(pat, tric_names[nTr], i)
                                 ] = dot_product(list(a), list(b)[0])
                    else:
                        dist_set[(pat, tric_names[nTr], i)
                                 ] = pearson_correlation(list(a), list(b)[0])
            nTr += 1
        pat += 1

    f_matrix = list()
    p = 0
    for bic_p in patients_tri:
        vals = filter(lambda x: x[0] == p, dist_set.keys())
        line = [dist_set[v] for v in vals]
        f_matrix.append(line)
        p += 1

    vals = filter(lambda x: x[0] == 0, dist_set.keys())
    cols = list(map(lambda z: str(z[1]) + "_" + str(z[2]), vals))

    return f_matrix, cols


def write_matrix(final_matrix, file_name, classes, cols, target, prog_r=None):
    fin = pd.DataFrame(data=final_matrix, columns=cols)
    fin[target] = classes
    fin.to_csv(file_name)


def write_matrix_static(static_feats_matrix, final_matrix, file_name, classes):
    cols = ["Gender", "BMI", "MND familiar history", 'Age at onset'] + \
        ["TRI"+str(n) for n in range(len(final_matrix[0]))]
    final_final_matrix = np.concatenate(
        (static_feats_matrix, final_matrix), axis=1)
    fin = pd.DataFrame(data=final_final_matrix, columns=cols)
    fin["Evolution"] = classes
    fin.to_csv(file_name)


######## AUXILIARES ##########

def pearson_correlation(numbers_x, numbers_y):

    for a, b in list(zip(numbers_x, numbers_y)):
        if a == math.nan or b == math.nan:
            return math.nan

    mean_x = sum(numbers_x) / len(numbers_x)
    mean_y = sum(numbers_y) / len(numbers_y)

    subtracted_mean_x = [i - mean_x for i in numbers_x]
    subtracted_mean_y = [i - mean_y for i in numbers_y]

    subtracted_mean_square_x = [(i - mean_x)**2 for i in numbers_x]
    subtracted_mean_square_y = [(i - mean_y)**2 for i in numbers_y]

    x_times_y = [a * b for a,
                 b in list(zip(subtracted_mean_x, subtracted_mean_y))]

    dem = (math.sqrt(sum(subtracted_mean_square_x))
           * math.sqrt(sum(subtracted_mean_square_y)))

    if dem == 0:
        return 0

    result = sum(x_times_y) / dem

    return result


def dot_product(lst1, lst2):

    a1, a2 = lst1
    b1, b2 = lst2

    if a1 == math.nan or a2 == math.nan or b1 == math.nan or b2 == math.nan:
        return math.nan

    inner_product = a1*b1 + a2*b2

    len1 = math.hypot(a1, a2)
    len2 = math.hypot(b1, b2)

    dem = len1*len2
    if dem == 0:
        return 0

    result = inner_product/dem

    return result


if __name__ == "__main__":
    ts = explore_triclusters("../outputs_3TPS/out1.txt")

    # for t in ts:
    #     print("Times:", t.times)
    #     print("Samples:", t.samples)
    #     print("Patients:", t.patients)
