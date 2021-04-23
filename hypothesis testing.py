from scipy import stats
from scipy.stats import ttest_rel, ttest_ind

from statsmodels.stats.multitest import multipletests
import numpy as np
from tabulate import tabulate


# Extract the file "lusc-rsem-fpkm-tcga_paired.txt" that contains the info of healthy genes into a list
list_of_paired = []
Paired = open("lusc-rsem-fpkm-tcga_paired.txt")
for line in Paired:
    # to remove any unwanted space at the begining or the end of the line
    stripped_line = line.strip()
    line_list = stripped_line.split()  # to put elements in a list
    list_of_paired.append(line_list)
Paired.close()
# print(len(list_of_paired))
# print(list_of_paired)

# Extract the file "lusc-rsem-fpkm-tcga-t_paired.txt" that contains the info of cancerous genes into a list
list_of_t_paired = []
Paired_t = open("lusc-rsem-fpkm-tcga-t_paired.txt")
for line in Paired_t:
    # to remove any unwanted space at the begining or the end of the line
    stripped_line = line.strip()
    line_list = stripped_line.split()  # to put elements in a list
    list_of_t_paired.append(line_list)
Paired_t.close()
# print(list_of_t_paired)

# lists of names of genes in healthy and cancerous cases
Healthy_names = []
Canser_names = []

# filtiration of genes
healthy_without_Zeros = []
cancer_without_Zeros = []

for a, b in zip(list_of_paired[1:], list_of_t_paired[1:]):
    Zeros = 0
    Zeroz = 0

    for x, y in zip(a, b):
        if x == '0.0':
            Zeros += 1
        if y == '0.0':
            Zeroz += 1
    if Zeros <= 25 and Zeroz <= 25:
        Healthy_names.append(a[0])
        Canser_names.append(b[0])
        healthy_without_Zeros.append(a[2:])
        cancer_without_Zeros.append(b[2:])
healthy_without_Zeros = [list(map(float, sublist))
                         for sublist in healthy_without_Zeros]
cancer_without_Zeros = [list(map(float, sublist))
                        for sublist in cancer_without_Zeros]
# print(cancer_without_Zeros)
# print(len(cancer_without_Zeros))
# print(Canser_names)

# get the p-values when Samples are paired.
p_val_pair = []
for x, y in zip(healthy_without_Zeros, cancer_without_Zeros):
    l = stats.ttest_rel(x, y).pvalue
    p_val_pair.append(l)


# get the p-values when Samples are independent.
p_val_ind = []
for x, y in zip(healthy_without_Zeros, cancer_without_Zeros):
    l = stats.ttest_ind(x, y).pvalue
    p_val_ind.append(l)



# #Apply the FDR multiple tests correction method on the paired samples.
# #list of tubles(reject:true for hypothesis that can be rejected for given alpha,
# #p-values corrected,corrected alpha for Sidak method,corrected alpha for Bonferroni method)
corrected_p_valpair_rej = multipletests(
    p_val_pair, alpha=0.05, method='fdr_bh')[0]
corrected_p_val_pair = multipletests(
    p_val_pair, alpha=0.05, method='fdr_bh')[1]
# print(corrected_p_val_pair)
# print(corrected_p_val_pair_rej)


# Apply the FDR multiple tests correction method on the independent samples.
corrected_p_val_ind_rej = multipletests(
    p_val_ind, alpha=0.05, method='fdr_bh')[0]
corrected_p_val_ind = multipletests(p_val_ind, alpha=0.05, method='fdr_bh')[1]
# print(corrected_p_val_ind_rej)


# #significance in the paired samples before FDR
# list of true for hypothesis that can be rejected for given alpha

# print(sig_p_val_pair_rej)
# print(len(sig_p_val_pair_rej))

# #get common significant genes for paired samples
sig_p_val_pair_com = []  # p-values of significant genes before fdr
sig_p_val_fdr_pair_com = []  # p-values of significant genes after fdr

# #get distinct significant genes for paired samples
sig_p_val_pair_dis = []
sig_p_val_fdr_pair_dis = []
for i in range(len(p_val_pair)):
    if p_val_pair[i] < 0.05 and corrected_p_val_pair[i] >= 0.05:
        sig_p_val_pair_dis.append(p_val_pair[i])
        sig_p_val_fdr_pair_dis.append(corrected_p_val_pair[i])
    elif corrected_p_val_pair[i] < 0.05:
        sig_p_val_pair_com.append(p_val_pair[i])
        sig_p_val_fdr_pair_com.append(corrected_p_val_pair[i])

# print(sig_p_val_pair_com)
print(len(sig_p_val_fdr_pair_com))

# print(sig_p_val_ind_com)
print(len(sig_p_val_fdr_pair_dis))

# #significance in the independent samples before FDR

# #get common significant genes for independent samples
sig_p_val_ind_com = []  # p-values of significant genes before fdr
sig_p_val_fdr_ind_com = []  # p-values of significant genes after fdr

# #get distinct significant genes for independent samples
sig_p_val_ind_dis = []
sig_p_val_fdr_ind_dis = []
for i in range(len(p_val_ind)):
    if p_val_ind[i] < 0.05 and corrected_p_val_ind[i] >= 0.05:
        sig_p_val_ind_dis.append(p_val_ind[i])
        sig_p_val_fdr_ind_dis.append(corrected_p_val_ind[i])
    elif corrected_p_val_ind[i] < 0.05:
        sig_p_val_ind_com.append(p_val_ind[i])
        sig_p_val_fdr_ind_com.append(corrected_p_val_ind[i])


# print(sig_p_val_ind_com)
print(len(sig_p_val_fdr_ind_com))
# print(sig_p_val_ind_com)
print(len(sig_p_val_fdr_ind_dis))

print(len(sig_p_val_fdr_ind_dis))
#table for independent case(common genes in first row vs distinct in the second row)
table1 = [sig_p_val_fdr_ind_com[0:6],sig_p_val_fdr_ind_dis[0:6]]
print(tabulate(table1))
#table for independent case(common genes in first row vs distinct in the second row)
table2 = [sig_p_val_fdr_pair_com[0:6],sig_p_val_fdr_pair_dis[0:6]]
print(tabulate(table2))