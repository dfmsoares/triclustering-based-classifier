# Learning Predictive Models with a Triclustering-based Classifier

This work proposes a new class of explainable predictive models from three-way data using discriminative triclusters. Triclustering searches are dynamically hyperparameterized to comprehensively find for informative triclusters (groups of individuals with a coherent pattern along a subset of variables and time points) and then use these temporal patterns as features within a state-of-the-art classifier with strict guarantees of interpretability. To this end, a temporally constrained triclustering algorithm, termed `TCtriCluster` algorithm, is devised to mine time-contiguous triclusters.

## Getting Started

These instructions will get you how to run the triclustering-based classifier. See the demo section for notes on how to run the demo project on your machine.

### Prerequisites

To run the triclustering-based classifier you need to have Python 3.4 or above installed as well as the following packages:

- [scikit-learn](https://scikit-learn.org/stable/install.html)
- [scipy](https://scipy.org/install.html)
- [numpy](https://numpy.org/install/)
- [pandas](https://pandas.pydata.org/getting_started.html)
- [sortedcontainers](http://www.grantjenks.com/docs/sortedcontainers/#quickstart)

### How To Run

First, you should perform triclustering in your data.

Run the following command to see the arguments list for `TCtriCluster`.

```
python3 TCtriCluster.py -h
```

Run the triclustering algorithm with the your defined input parameters and save the result for an output file (`.txt`):

```
python3 TCtriCluster.py -f <input_file> -sT <min_t> -sS <min_s> -sG <min_g> -w <win_ratio> -o <opc> [-mv <mv_threshold>] > <output_file>.txt
```

Next, with the produced triclusters, you compute the similarity matrices:

```
python3 compute_similar_mats_tri.py <datafile> <triclusters_file> <matrices_folder> <target_column> <nr_time_points> <categorical_features> <continuos_features>
```

Finally, the classifier uses the similarity matrices as the learning examples. This code performs `n` x `k-fold` Stratified CV to evaluate the performance of the classifier.

```
python3 compute_predictions.py <matrices_folder> <target_column> <nr_time_points> <output_csv> <k_splits> <n_repeats>
```

## Demo Example

We provide a demo example in [`demo`](/demo) folder.

First we parsed the [d1.csv](/demo/d1.csv) file to achieve the required formatting by the triclustering algorithm. To do this, use the following command:

```
# Usage python3 src/write_tab_file.py <csv_file> <output_tabfile> <target_var> <n_timepoints>

$ python3 src/write_tab_file.py demo/d1.csv demo/tab_file.tab Class 3
```

Next we performed triclustering with TCtriCluster algorithm:

```
$ python3 src/TCtriCluster.py -f demo/tab_file.tab -sT 2 -sS 2 -sG 2 -w 0.1 -o 1 > demo/triclusters_d1.txt
```

Next, with triclusters and the original dataset we computed the similiarity matrices

```
$ python3 src/compute_similar_mats_tri.py demo/d1.csv demo/triclusters_d1.txt demo/sim_matrices Class 3 [S1, S2, S3, S4, S5, S6] []
```

Finnally with the matrices we run the classifier, evaluating the results with repeated stratified k-fold CV:

```
$ python3 src/compute_predictions.py demo/sim_matrices Class 3 results.csv 2 2
```

## Citing the Paper 📑

If you use the triclustering-based classifier in your research, please cite our paper:

*Soares, D. F., Henriques, R., Gromicho, M., de Carvalho, M., & C Madeira, S. (2023) Triclustering-based classification of longitudinal data for prognostic prediction: targeting relevant clinical endpoints in amyotrophic lateral sclerosis. Scientific Reports, 6182* https://doi.org/10.1038/s41598-023-33223-x

