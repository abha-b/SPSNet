# SPSNet
Code of  SPSNet: Subpopulation-sensitive network-based analysis of heterogeneous gene expression data. BMC Systems Biology, 12(Suppl 2):28, April 2018

*************************************************************************************************
```
Program description: 
	This program implements SPSNet, PFSnet as described in the paper -- 
	Abha Belorkar, Rajanikanth Vadigepalli, and Limsoon Wong. 
	"SPSNet: subpopulation-sensitive network-based analysis of 
	heterogeneous gene expression data" 
	BMC Systems Biology. 2018.
Date created: 10/04/2015
Last modified: 23/02/2018
Author: Abha Belorkar
Contact: abha.b.nus@gmail.com, wongls@comp.nus.edu.sg, Rajanikanth.Vadigepalli@jefferson.edu
Affiliated organizations: 
	School of Computing, National University of Singapore, Singapore
	Department of Pathology, Thomas Jefferson University, Philadelphia
```
	
*************************************************************************************************

```
REQUIREMENTS:
	Python 2.7
	packages required: numpy, scipy, networkx, rpy2, sklearn, seaborn, mpl_toolkits

CALL:
	python spsnet.py
	
ARGS:
	1. --control=<control_expr_file> OR -c <control_expr_file> (compulsory argument)
			Absolute/relative path of the file containing expression matrix for control group samples; tab separated text file; matrix: genes x samples (first column: gene ids/names); first line skipped as header.
	2. --test=<test_expr_file> OR -t <test_expr_file> (compulsory argument)
			Absolute/relative path of the file containing expression matrix for test group samples; tab separated text file; matrix: genes x samples (first column: gene ids/names); first line skipped as header.
	3. --pathway=<pathway_file> OR -p <pathway_file> (compulsory argument)
			Tab separated txt file; first line skipped as header; each edge as a row; 3 columns - pathway name<tab>gene_1<tab>gene_2.
	4. --pathway=<pathway_file> OR -p <pathway_file> (compulsory argument)
			Tab separated txt file; first line skipped as header; each edge as a row; 3 columns - pathway name<tab>gene_1<tab>gene_2.
	5. --theta1=<theta1_value> OR -h <theta1_value> (optional; default is 0.95)
			value between 0 to 1; denotes the quantile threshold above which patient gives full vote to a gene.
	6. --theta2=<theta2_value> OR -l <theta2_value> (optional; default is 0.85)
			value between 0 to 1; denotes the quantile threshold below which patient gives zero vote to a gene.
	7. --beta=<beta_value> OR -b <beta_value> (optional; default is 0.5)
			value between 0 to 1; denotes the minimum average vote for which the gene is considered highly expressed.\n\n"
	8. --option=<option_var> OR -o <option_value> (optional; default is SPS)
			SPS or PFS; denotes the method of differentially expression analysis.
	9. --typesize=<typesize_value> OR -k <typesize_value> (optional; default is 10)
			value greater than 5; denotes the number of samples selected to represent potential subpopulations.
	10. --symmetry=<symmetry_var> OR -s <symmetry_var> (optional; default is 'n')
			y or n; denotes whether heterogeneity in control group should also be accounted for 
			Value y is recommended for analyzing datasets with multiple batches.
	
EXAMPLE CALLS:
		python spsnet.py -c datasets/ALL/training/control.txt -t datasets/ALL/training/subtype1_30_subtype2_29.txt -p pathways/pathwayAPI_human_pathways_entrez_id.txt -o SPS -i datasets/ALL/training/S1_30_S2_29_sample_info.txt
		python spsnet.py -c datasets/HCC/control_expr_entrez.txt -t datasets/HCC/hcc_expr_entrez.txt -p pathways/pathwayAPI_human_pathways_entrez_id.txt -o SPS -i datasets/HCC/sample_info.txt -s y
		python spsnet.py -c datasets/rat_toxicogenomics/PPARA_MoA/GSE55347_ge_control.txt -t datasets/rat_toxicogenomics/PPARA_MoA/GSE55347_ge_test.txt -p pathways/KEGG_rat_pathways_gene_symbols.txt -o SPS -i datasets/rat_toxicogenomics/PPARA_MoA/sample_info.txt
		
EXPECTED OUTPUT:
	A folder named 'results/PFS/' or 'results/SPS/' is created at the location where the script is run (depending on the -o argument, default: SPS). This contains:
	a. sps_b_0.5: A list of subnetworks reported significant by SPSNet/PFSNet, their p-values, and the genes forming the subnetworks
	b. pc_loadings.txt: PC loadings (first 3 PCs) corresponding to each significant subnetwork
	c. sps_pc12.png, sps_pc23.png, sps_pc13.png, sps_pc123.png: PCA scatter plots generated using SPS/PFS scores of significant subnetworks
```	
*************************************************************************************************
