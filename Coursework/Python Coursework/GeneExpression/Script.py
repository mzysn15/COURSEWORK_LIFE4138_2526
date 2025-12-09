############################################################
############################################################
 #🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬#
 #🧬                                                    🧬#
 #🧬         This is my Python Script sumbmission       🧬#
 #🧬       for the Gene Expression Dataset, Set 12.     🧬#
 #🧬                                                    🧬#
 #🧬                                                    🧬#
 #🧬              Author: Shahwar Nadeem                🧬#
 #🧬                                                    🧬#
 #🧬                                                    🧬#
 #🧬   FYI: Code below is for 1/2 datasets, not both,   🧬#
 #🧬    but the same code applies for both datasets,    🧬#
 #🧬       and has been used on both data sets!         🧬#
 #🧬                        🫶                          🧬#
 #🧬                                                    🧬#
 #🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬🧬#
############################################################
############################################################

#Import ALL Libraries
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import math
import Bio


#Read the .tsv files into python and display them
DF_A = pd.read_csv("A_vs_E.deseq2.results.tsv", sep = '\t')
print(DF_A)

#DF_A - Basic Summary of P-values and Log Fold Changes
DF_A[['log2FoldChange', 'pvalue']].describe()

#DF_A - Identify signficant Up/Downregulation in genes. If fold change above 0 + P-adj below 0.05 then return "Up". If fold change below 0 + P-adj above 0.05 then return "Down".
conditions = [
    (DF_A['log2FoldChange'] > 0) & (DF_A['padj'] < 0.05),
    (DF_A['log2FoldChange'] < 0) & (DF_A['padj'] < 0.05)]
choices = ['Up', 'Down']

DF_A['Up/Down_Reg'] = np.select(conditions, choices, default='NotSign')

#Significantly Upregulated and Downregulated genes with their corresponding fold changes, p-values, and adjusted p-values.
print(DF_A[['gene_id', 'log2FoldChange', 'pvalue', 'padj', 'Up/Down_Reg']].head())

#Remove all Rows with null values
DF_A = DF_A.dropna()

#Count the number of Up/Down Regulations and Not Significant in DF_A, Significant Gene List
DF_A.groupby('Up/Down_Reg').count()
##Significantly Upregulated and Downregulated genes with their corresponding fold changes, p-values, and adjusted p-values.
DF_A[['gene_id', 'log2FoldChange', 'pvalue', 'padj', 'Up/Down_Reg']]
##Significantly Upregulated genes will display only
UpReg = DF_A[DF_A['Up/Down_Reg'] == 'Up'][['gene_id', 'log2FoldChange', 'pvalue', 'padj', 'Up/Down_Reg']]
UpReg
##Significantly Downregulated genes will display only
DownReg = DF_A[DF_A['Up/Down_Reg'] == 'Down'][['gene_id', 'log2FoldChange', 'pvalue', 'padj', 'Up/Down_Reg']]
DownReg

#Summary of p-values and log fold changes across all genes for each comparison.
x = (DF_A['log2FoldChange'])
y = (DF_A['pvalue'])
colors = np.random.randint(100, size = len(DF_A))
plt.scatter(x, y, c = colors, cmap = 'Blues')
plt.colorbar()
plt.xlabel("Log 2 Fold Change", fontsize=10)
plt.ylabel("P Value", fontsize=10)
plt.title("Summary of P values and Log Fold Changes across all Genes")
plt.show()

#DF_A - Volcano plots to visualize the significance and magnitude of changes in gene expression for each comparison. I fear I may have gone the long way around this...but this was only tutorial where I fully understood what each code was doing (maybe not why)!
##Calculate log10 Value
DF_A['NegLog_Pval'] = -np.log10(DF_A['pvalue'])
plt.figure(figsize = (12,8))

##Scatter plot of Upregulated Genes
plt.scatter(DF_A['log2FoldChange'][(DF_A['padj'] < 0.05) & (DF_A['log2FoldChange'] >= 2.5) & (DF_A['NegLog_Pval'] >= 5)],
            DF_A['NegLog_Pval'][(DF_A['padj'] < 0.05) & (DF_A['log2FoldChange'] >= 2.5) & (DF_A['NegLog_Pval'] >= 5)],
c = 'red', label = 'Upregulated')
##Scatter Plot of Downregulated Genes
DF_A['NegLog_Pval'] = -np.log(DF_A['pvalue'])
plt.scatter(DF_A['log2FoldChange'][(DF_A['padj'] < 0.05) & (DF_A['log2FoldChange'] <= -2.5) & (DF_A['NegLog_Pval'] >= 5)],
            DF_A['NegLog_Pval'][(DF_A['padj'] < 0.05) & (DF_A['log2FoldChange'] <= -2.5) & (DF_A['NegLog_Pval'] >= 5)],
c = 'blue', label = 'Downregulated')

#Volcano Plot of Up/Downregulated Genes in DF_A
for i, gene in enumerate (DF_A['gene_id']):
    plt.annotate(gene, (DF_A['log2FoldChange'][i], DF_A['NegLog_Pval'][i]))
plt.title("Volcano Plot", fontsize = 20)
plt.xlabel("Log2 Fold Change", fontsize = 10)
plt.ylabel("-Log10 P-Value", fontsize = 10)
plt.legend(loc = 'upper right', fontsize = 8)
plt.axvline(5, ls = '--', color = 'black')
plt.axvline(2.5, ls = '--', color = 'black')
plt.axvline(-2.5, ls = '--', color = 'black')
plt.show()

#DF_A - MA Plot to display the relationship between the Log Fold change and Mean Expression for all genes, and significant genes highlighted
DF_A['Log_BaseMean'] = np.log10(DF_A['baseMean'])
plt.figure(figsize = (12, 8))
plt.scatter(DF_A['Log_BaseMean'], DF_A['log2FoldChange'],
            color = 'grey', s = 20, label = 'Gene')
plt.scatter(DF_A.loc[DF_A['padj'] < 0.05, 'Log_BaseMean'],
            DF_A.loc[DF_A['padj'] < 0.05, 'log2FoldChange'],
            color = 'black', s = 20, label = 'Signifcant Genes')
plt.title("Relationship between the Log Fold change and Mean Expression", fontsize = 20)
plt.xlabel("log10 Base Mean", fontsize = 10)
plt.ylabel("Log2 Fold Change", fontsize = 10)
plt.axhline(0, ls = '--', color = 'black')
plt.legend(loc = 'upper right', fontsize = 8)
plt.show()

#DF_A - Histogram of P-value Distribution
plt.hist(DF_A['pvalue'], bins = 25, color = 'white', edgecolor = 'black')
plt.title("Distribution of P-values")
plt.xlabel("P-values")
plt.ylabel("Frequency")
plt.show()

#Heatmap for top differentially expressed genes across conditions-3 Heatmaps will be produced. HM1 will display top differentialy expressed genes for Condition A vs E. HM2 will display top differentially expressed genes for Condition I vs E. HM3 will display the top differentially expressed genes common to both conditions but their expression is different across both conditions. 

#HEATMAP 1-Top Differentialy Expressed Genes in A vs E Deseq2 Results
DF_A1v = DF_A.sort_values(by=["log2FoldChange", "padj"],
                             ascending=[False, True])
DF_A1v = DF_A1v.rename(columns={'log2FoldChange': 'log2fc_A'})
DF_A1v = DF_A1v.head(15)
sns.heatmap(DF_A1v[['log2fc_A']].set_index(DF_A1v['gene_id']), annot = True, fmt=".1f", cmap= 'viridis')
plt.title("The top differentially expressed genes", fontsize = 15, weight = "bold")
plt.ylabel("Gene ID", fontsize = 10)
plt.xlabel("Log 2 Fold Change in A vs E Dataset")
plt.show()

#HEATMAP 2-Top Differentialy Expressed Genes in I vs E Deseq2 Results
DF_I = pd.read_csv("I_vs_E.deseq2.results.tsv", sep = '\t')

regulation = [
    (DF_I['log2FoldChange'] > 0) & (DF_I['padj'] < 0.05),
    (DF_I['log2FoldChange'] < 0) & (DF_I['padj'] < 0.05)]
choices = ['Up', 'Down']

DF_A['Up/Down_Reg'] = np.select(regulation, choices, default='NotSign')
print(DF_I)

DF_I1v = DF_I.sort_values(by=["log2FoldChange", "padj"],
                             ascending=[False, True])
DF_I1v = DF_I1v.rename(columns={'log2FoldChange': 'log2fc_I'})
DF_I1v = DF_I1v.head(15)
sns.heatmap(DF_I1v[['log2fc_I']].set_index(DF_I1v['gene_id']), annot = True, fmt=".1f", cmap= 'viridis')
plt.title("The top differentially expressed genes", fontsize = 15, weight = "bold")
plt.ylabel("Gene ID", fontsize = 10)
plt.xlabel("Log 2 Fold Change in I vs E Dataset")
plt.show()

#HEATMAP 3-Top Differentialy Expressed Common Genes across the Conditions
DFA3v = DF_A.sort_values(by=["log2FoldChange", "padj"],
                             ascending=[False, True])
DFA3v = DFA3v.head(100)

DFI3v = DF_I.sort_values(by=["log2FoldChange", "padj"],
                             ascending=[False, True])
DFI3v = DFI3v.head(100)

DFA3v = DFA3v.rename(columns={'log2FoldChange': 'log2fc_A'})
DFI3v = DFI3v.rename(columns={'log2FoldChange': 'log2fc_I'})

MERGED4v = DFA3v.merge(DFI3v, on = 'gene_id')
sns.heatmap(MERGED4v[['log2fc_I', 'log2fc_A']].set_index(MERGED4v['gene_id']), annot = True, fmt=".1f", cmap= 'viridis')
plt.title("The top differentially expressed genes", fontsize = 15, weight = "bold")
plt.ylabel("Gene ID", fontsize = 10)
plt.show()

#Additional Analyses-This is a Clustered Heatmap
sns.clustermap(MERGED4v[['log2fc_I', 'log2fc_A']].set_index(MERGED4v['gene_id']), 
            annot=True, figsize=(7,6), dendrogram_ratio=(.6, .05), cbar_pos=(1, .15, .04, .7))
plt.show()




#I am concerned I have not used def or _init_ at all in this code, but I do not know why I have to use those codes, I cannot imagine using def to define my graphs when i can use '#...' to explain what is happening. Furthermore, the sites which have supplemented my learning/coursework are: GeeksforGreeks, W3 Schools, @Bioinformatics-IBE, @MrBioinformatiX, Seaborn (+).







