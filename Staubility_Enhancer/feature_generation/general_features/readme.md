These functions generate various features related to the properties of yeast's genes.
This readme summarizes the features calculated within this directory, as well as the inputs and outputs of all functions.

Hic:
Hic is a matrix of 'smoothed Hi-C graph distance between S. cerevisiae genes'.
The function calculates for each gene:
1.	Rank - number of connections
2.	Mean value of distances
3.	Median value of distances 
4.	Min value of distances
5.	Mean mRNA levels of each gene and its connected genes

inputs - HiC.mat , hic_data.csv 

output - HiC.csv

CoExp:
CoExp.mat calculates a correlation score between each pair of genes according to their mRNA levels, using different databases.
This function calculates the weighted rank for each gene according to the CoExp matrix.

files needed - CoExp.mat 

output - Weighted_Rank.csv

lfe:
This function generates the mRNA folding energy features.
These features are based on the following tables, which contain lfe values in different windows.
The first table is - native local folding energy (empirical results).
The second table is - shuffled lfe (based on a random model).
The third table is - dlfe, the difference between the results in table 1 and table 2.
In addition, we calculated the mean value per gene for each table.

inputs - export_taxid_559292_begin_native.csv, export_taxid_559292_begin_shuffled.csv, export_taxid_559292_begin_dlfe.csv.
output - dlfe.csv, lfe_native.csv, lfe_shuffled.csv, mean_lfe.csv

sORF:
sORFs are shifted Open Reading Frames, beginning at alternative ATGs downstream in the ORF from the main START codon, and terminating with a stop codon.
The features we calculated using the gene's sequence:
1.	Number of shifted ORFs
2.	Number of shifted ORFs 30 codons
3.	Number of shifted ORFs 200 codons
4.	Max length
5.	Max length 30 codons
6.	Max length 200 codons
7.	Mean length
8.	Mean length 30 codons
9.	Mean length 200 codons

inputs - full_table.csv (contains the sequence of each gene)
output - sORF.csv, sORF_Normalized.csv

ATG_features:
We used the table 'ATG_Tamir' which contains ATG context score for 5861 genes.
The score was calculated for the main ATG, for each alternative ATG in the 150 nucleotides of the untranslated region on the 5’ side (5UTR) and for each alternative ATG in the 150 nucleotides of the open reading frame (ORF).
For each alternative ATG, an absolute score and a relative score (absolute score divided by the main ATG score) were  calculated.
The features calculated:
i.	The number of alternative ATG in the following windows (default value=0):
1.	ORF
2.	5UTR
3.	Total
4.	30 codons in the ORF
5.	200 codons in the ORF
Alternative ATGs can affect translation and impair it.
ii.	Main ATG context score.
iii.	Mean and maximum values for the scores in the following windows (default value=10 orders of magnitude less than the worst score):
1.	30 codons in the ORF
These features were calculated for both absolute and relative scores.

inputs - full_table.csv, ATG_Tamir.xlsx
output - ATG_withnan1.csv, ATG_withnan.csv

In addition we used the 'Gravy_aliphatic' and 'fop_instability' tables, containing for each yeast's gene the following features:

Aliphatic index- The aliphatic index of a protein is defined as the relative volume occupied by aliphatic side chains (alanine, valine, isoleucine, and leucine). Aliphatic index plays a role in protein thermal stability. Proteins with a high Aliphatic index are more thermally stable. Aliphatic amino acids also are hydrophobic in nature.

Gravy score- The GRAVY number of a protein is a measure of its hydrophobicity or hydrophilicity. The two measures are combined in a hydropathy scale or hydropathy index. The more positive the value, the more hydrophobic the amino acids in this region. These scales are commonly used to predict the transmembrane alpha helices of membrane proteins. 

Instability index- The Instability index is a measure of proteins, used to determine whether it will be stable in a test tube. If the index is less than 40, then it probably is. .
Fop- This index measures the frequency of optimal codons (Fop) in a gene. It is a species-specific measure of bias towards certain codons that appear to be translationally optimal in a species.



mail info: glicksteinbar@gmail.com
