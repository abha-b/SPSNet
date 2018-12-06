#! /usr/bin/python

'''
Description: 
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
'''

import os
import sys, getopt
import numpy as np
from scipy.stats import rankdata
import networkx as nx
from time import time
import rpy2.robjects as R
from numpy.random import choice
from math import sqrt
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

P_VAL_CUT_OFF = 0.05

# prints program usage instructions - command and expected arguments
def usage():
	
	print "\nCALL:\n"
	print "python spsnet.py\n\n"
	
	print "ARGS:\n"
	print "1. --control=<control_expr_file> OR -c <control_expr_file> (compulsory argument)\n"
	print "Absolute/relative path of the file containing expression matrix for control group samples; tab separated text file; matrix: genes x samples (first column: gene ids/names); first line skipped as header.\n\n"
	print "2. --test=<test_expr_file> OR -t <test_expr_file> (compulsory argument)\n"
	print "Absolute/relative path of the file containing expression matrix for test group samples; tab separated text file; matrix: genes x samples (first column: gene ids/names); first line skipped as header.\n\n"
	print "3. --pathway=<pathway_file> OR -p <pathway_file> (compulsory argument)\n"
	print "Tab separated txt file; first line skipped as header; each pathway edge as a row; 3 columns - pathway name<tab>gene_1<tab>gene_2.\n\n"
	print "4. --info=<sample_info_file> OR -i <sample_info_file> (compulsory argument)\n"
	print "Tab separated txt file; first line skipped as header; each sample info as a row; 2 columns - sample_name<tab>class.\n\n"
	print "5. --theta1=<theta1_value> OR -h <theta1_value> (optional; default is 0.95)\n"
	print "value between 0 to 1; denotes the quantile threshold above which patient gives full vote to a gene.\n\n"
	print "6. --theta2=<theta2_value> OR -l <theta2_value> (optional; default is 0.85)\n"
	print "value between 0 to 1; denotes the quantile threshold below which patient gives zero vote to a gene.\n\n"
	print "7. --beta=<beta_value> OR -b <beta_value> (optional; default is 0.5)\n"
	print "value between 0 to 1; denotes the minimum average vote for which the gene is considered highly expressed.\n\n"
	print "8. --option=<option_var> OR -o <option_value> (optional; default is SPS)\n"
	print "SPS or PFS; denotes the method of differentially expression analysis.\n\n"
	print "9. --typesize=<typesize_value> OR -k <typesize_value> (optional; default is 10)\n"
	print "value greater than 5; denotes the number of samples selected to represent potential subpopulations.\n\n"
	print "10. --symmetry=<symmetry_var> OR -s <symmetry_var> (optional; default is 'n')\n"
	print "y or n; denotes whether heterogeneity in control group should also be accounted for.\n\n"
	
	sys.exit(2)

# accepts arguments, performs necessary checks, assigns defaults
def handle_args (argv):

	# initialization
	fe_c_name = ""
	fe_t_name = ""
	fpw_name = ""
	fs_i_name = ""
	beta = 0.5
	theta_1 = 0.95
	theta_2 = 0.85
	k = 5
	option = 'SPS'
	symmetry = 'N'
	sel_indices = 0

	# processing
	try:
		opts, args = getopt.getopt(argv, "c:t:p:i:b:h:l:n:o:k:s:", \
			["control=","test=","pathway=","info=","beta=","theta1=","theta2=","option=",\
			"typesize=","symmetry="])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-c","--control"):
			fe_c_name = arg
		elif opt in ("-t","--test"):
			fe_t_name = arg
		elif opt in ("-p","--pathway"):
			fpw_name = arg
		elif opt in ("-i","--info"):
			fs_i_name = arg
		elif opt in ("-b","--beta"):
			beta = float(arg)
		elif opt in ("-h","--theta1"):              
			theta_1 = float(arg)
		elif opt in ("-l","--theta2"):
			theta_2 = float(arg)
		elif opt in ("-o","--option"):
			option = arg.upper()
		elif opt in ("-k","--typesize"):
			k = int(arg)
		elif opt in ("-s","--symmetry"):
			symmetry = arg.upper()[0]
		else:
			print "Unrecognized argument: ", opt, "\n\n"
			usage()

	# handling empty arguments
	if not (fe_c_name):
		print "\nEnter the path of control group gene expression file (with -c).\n"
		usage()
	if not (fe_t_name):
		print "\nEnter the path of test group gene expression file (with -t).\n"
		usage()
	if not (fpw_name):
		print "\nEnter the path of pathway file (with -p).\n"
		usage()
	if not (fs_i_name):
		print "\nEnter the path of sample info file (with -i).\n"
		usage()
	if not (beta < 1 and beta > 0):
		print "\nEnter a beta value between 0 and 1.\n"
		usage()
	if not (theta_1 < 1 and theta_1 > 0):
		print "\nEnter a theta1 value between 0 and 1.\n"
		usage()
	if not (theta_2 < 1 and theta_2 > 0):
		print "\nEnter a theta2 value between 0 and 1.\n"
		usage()	
	if not (option == 'PFS' or option == 'SPS'):
		print '\nEnter an appropriate option: PFS or SPS.\n'
		usage()
	if not (k >= 5):
		print '\nEnter a k value greater than 5.\n'
		usage()
	if not (symmetry == 'Y' or symmetry == 'N'):
		print '\nEnter y or n.\n'
		usage()
		
	return fe_c_name, fe_t_name, fs_i_name, fpw_name, theta_1, theta_2, beta, k, option, symmetry

# returns gene expression matrix
def load_file(fi_name):

	labels, data = [], []
	fi = open(fi_name, 'r')
	snames = np.array(fi.readline().strip().split('\t')[1:])
	skipped_gcount = 0

	# file --> 2-d array (genes x samples)
	# skip lines (genes) with missing values
	for line in fi:
		ln = line.strip().split('\t')
		try:
			row = map(float, ln[1:])
			data.append(row)
			labels.append(ln[0].upper())
		except:
			# 'try' fails when missing values (blanks) cannot convert to float
			skipped_gcount += 1

	return np.array(data), labels, snames, skipped_gcount

# returns a graph per pathway
def load_pathways(fpw_name, genes):
	
	pw2G = {}

	# populate dict pw2G: key (pathway name), value (gene graph)
	fpw = open(fpw_name, 'rb')
	header = fpw.readline()
	for line in fpw:
		row = line.strip().split('\t')
		if len(row) == 3:
			pw, node_1, node_2 = row
			if node_1 in genes and node_2 in genes:
				if pw not in pw2G:
					pw2G[pw] = nx.Graph()
				pw2G[pw].add_edge(node_1.upper(), node_2.upper())
	fpw.close()

	return pw2G

# computes weights of genes in patients
def compute_weights(X, genes, theta_1, theta_2):

	num_genes, num_samples = np.shape(X)

	# function calculating vote to a gene from a patient, given quantiles in the group
	def get_vote(r, q1, q2):
		if r >= q1:
			return 1
		if r >= q2:
			return float(r-q2)/(q1-q2)
		else:
			return 0

	# for each sample, calculate gene expression ranks: output matrix -- samples x genes
	'''
	Note on default behavior of function scipy.stats.rankdata: 
	Ranks of tied elements are assigned to be the average of the ranks that would have been assigned 
	to all the tied values had they been unequal. (This is also the case in R by default.)
	'''
	ranks = [rankdata(X[:,j])/float(num_genes) for j in range(num_samples)]

	# for each sample, calculate theta_1 and theta_2 quantiles for gene expression
	q1 = np.percentile(ranks, theta_1*100, axis=1)
	q2 = np.percentile(ranks, theta_2*100, axis=1)

	# 'weights' is a dict; gene (key) --> fuzzy vote in all samples (value: list)
	weights = {genes[i]:[get_vote(ranks[j][i], q1[j], q2[j]) for j in range(num_samples)] for i in range(num_genes)}

	return weights
	
# break graph into subnets with immediate neighbors
def break_into_subnets(pw2G):

	subnets = {}
	for pw in pw2G:
		pw_nodes = pw2G[pw].nodes()
		for node in pw_nodes:
			s_nodes = pw2G[pw].neighbors(node) + [node]
			if len(s_nodes) >= 5:
				subnets[node + '__' + pw] = s_nodes

	return subnets
	
def scatter_2d(root, fmatrix_t, labels, cls_colors, i, j, var, option, sc):

	pc_percent_var = {pi: format(var[pi], '.2f') for pi in range(np.shape(fmatrix_t)[1])}

	# scatter plot
	fig = plt.figure(figsize=(12, 6))
	ax = fig.add_subplot(111)
	for class_name, col in zip(set(labels), cls_colors):
		#print class_name
		x, y = fmatrix_t[labels == class_name, i], fmatrix_t[labels == class_name, j]
		p = ax.scatter(x, y, color=col, label=class_name, alpha=1, s=150)
		'''
		texts = []
		for (xi, yi, li) in zip(x, y, range(len(x))):
			texts.append(ax.text(xi, yi, li))
		adjust_text(texts, force_text=0.1, force_points=0.7, arrowprops=dict(arrowstyle="-|>", color='r', alpha=0.5))
		'''
	ax.set_xlabel('PC ' + str(i+1) + ' - ' + pc_percent_var[i] + '% variance', fontsize=14, labelpad=12)
	ax.set_ylabel('PC ' + str(j+1) + ' - ' + pc_percent_var[j] + '% variance', fontsize=14, labelpad=12)
	#ax.set_ylim(-6.5, 15)
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
	
	handles, labels = ax.get_legend_handles_labels()
	# sort both labels and handles by labels
	labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
	
	leg = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints = 1, fontsize=11, ncol=1, frameon=True)
	leg.get_frame().set_facecolor('1.0')
	#plt.title('PCA plot of subtypes')
	plt.tick_params(axis='both', which='major', labelsize=13, pad=10)
	plt.tick_params(axis='both', which='minor', labelsize=13, pad=10)
	plt.title('Silhouette score: '+str(sc)[:5] + '\n', fontsize=14)
	fig.savefig(root + '/' + option.lower() + '_pc' + str(i+1) + str(j+1) + '.png')
	plt.close()

	return
	
def scatter_3d(root, fmatrix_t, labels, cls_colors, i, j, k, var, option, sc):

	pc_percent_var = {pi: format(var[pi], '.2f') for pi in range(np.shape(fmatrix_t)[1])}

	# scatter plot
	fig = plt.figure(figsize=(25,15))
	ax = fig.add_subplot(111, projection='3d')
	for class_name, col in zip(set(labels), cls_colors):
		#print class_name
		x, y, z = fmatrix_t[labels == class_name, i], fmatrix_t[labels == class_name, j], fmatrix_t[labels == class_name, k]
		p = ax.scatter(x, y, z, color=col, label=class_name, alpha=1, s=300)
		
	ax.set_xlabel('PC ' + str(i+1) + ' - ' + pc_percent_var[i] + '% variance', fontsize=14, labelpad=12)
	ax.set_ylabel('PC ' + str(j+1) + ' - ' + pc_percent_var[j] + '% variance', fontsize=14, labelpad=12)
	ax.set_zlabel('PC ' + str(k+1) + ' - ' + pc_percent_var[k] + '% variance', fontsize=14, labelpad=12)
	#ax.set_ylim(-6.5, 15)
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	
	handles, labels = ax.get_legend_handles_labels()
	# sort both labels and handles by labels
	labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
	
	leg = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints = 1, fontsize=11, ncol=1, frameon=True)
	leg.get_frame().set_facecolor('1.0')
	#plt.title('PCA plot of subtypes')
	plt.tick_params(axis='both', which='major', labelsize=13, pad=10)
	plt.tick_params(axis='both', which='minor', labelsize=13, pad=10)
	plt.title('Silhouette score: '+str(sc)[:5] + '\n', fontsize=14)
	fig.savefig(root + '/' + option.lower() + '_pc' +str(i+1)+ str(j+1) + str(k+1) + '.png')
	plt.close()

	return
	
# PCA scatter plots
def plot_PCA(d_path, f_matrix, option, classes_c, classes_t, features):

	if not os.path.isdir(d_path):
		os.makedirs(d_path)
	if not os.path.isdir(d_path + '/' + option):
		os.makedirs(d_path + '/' + option)
	root = d_path + '/' + option

	sns.set_style('white')

	fmatrix = np.transpose(f_matrix)
	nc = 30 if 30 < len(features) else len(features)
	pca = PCA(n_components=nc)
	
	print np.shape(fmatrix)
	fmatrix_t = pca.fit_transform(fmatrix)
	print np.shape(fmatrix_t)
	var = pca.explained_variance_ratio_
	labels = np.array(list(classes_c) + list(classes_t))
	#samples = np.array(list(snames_c) + list(snames_t))
	#print labels
	
	comp = pca.components_
	#print np.shape(comp)
	si = [i for i in range(len(features)) if np.any(comp[:, i])]
	sel_subnets = np.array(features)[si]
	
	num_pc, num_feat = np.shape(pca.components_)
	fpc = open(root + '/pc_loadings.txt', 'wb')
	fpc.write('feature\t' + '\t'.join(['pc_' + str(x+1) + '_loading' for x in range(num_pc)]) + '\n')
	for i in range(num_feat):
		fpc.write(features[i] + '\t' + '\t'.join([str(y) for y in pca.components_[:, i]]) + '\n')
	fpc.close()
	
	#print labels
	sc_123 = silhouette_score(fmatrix_t[:, :3], labels, metric='euclidean')
	sc_12 = silhouette_score(fmatrix_t[:, :2], labels, metric='euclidean')
	sc_23 = silhouette_score(fmatrix_t[:, 1:3], labels, metric='euclidean')
	sc_13 = silhouette_score(fmatrix_t[:, [0,2]], labels, metric='euclidean')
	print 'silhouette score (PC 1, 2):', sc_12
	print 'silhouette score (PC 2, 3):', sc_23
	print 'silhouette score (PC 2, 3):', sc_13
	print 'silhouette score (PC 1, 2, 3):', sc_123
	
	#'''
	red_CB = "#b30000"
	orange_CB = "#ff7f00"
	blue_CB = "#377eb8"
	green_CB = "#006d2c"
	grey_CB = "#999999"
	yellow_CB = "#ffff33"
	brown_CB = "#a65628"
	pink_CB = "#f781bf"
	purple_CB = "#984ea3"
	l_blue_CB = '#a6cee3'
	l_red_CB = '#fb9a99'
	l_green_CB = '#b2df8a'
	l_orange_CB = '#fdbf6f'
	maroon_CB = '#800000'
	olive_CB = '#808000'
	aqua_CB = '#00FFFF'
	lime_CB = '#00FF00'
	teal_CB = '#008080'
	navy_CB = '#000080'
	#'''
	
	#cls_colors = ['green', 'mediumaquamarine', 'red', 'violet']
	cls_colors = [green_CB, red_CB, pink_CB, orange_CB, blue_CB, purple_CB, grey_CB, brown_CB, olive_CB, aqua_CB, lime_CB]
	
	scatter_2d(root, fmatrix_t, labels, cls_colors, 0, 1, var, option, sc_12)
	scatter_2d(root, fmatrix_t, labels, cls_colors, 1, 2, var, option, sc_23)
	scatter_2d(root, fmatrix_t, labels, cls_colors, 0, 2, var, option, sc_13)
	scatter_3d(root, fmatrix_t, labels, cls_colors, 0, 1, 2, var, option, sc_123)
	
	return sel_subnets

# this function computes scores for each patient-subnetwork pair in a given class
def get_subnet_scores (subnets, weights_1, weights_2, num_samples_1, num_samples_2, beta, k, option, classes_c, classes_t, symm):

	# scores is a dictionary which stores the scores (values) of all patients for each subnet (key)
	scores_1 = {s_name:[] for s_name in subnets}
	scores_2 = {s_name:[] for s_name in subnets}
	subnet2pval = {}
	sgnf_subnets = []
		
	# for each patient:
	# score 1: sum over --> (avg fuzzy vote of gene in type 1) x (fuzzy vote to gene from patient)
	# score 2: sum over --> (avg fuzzy vote of gene in type 2) x (fuzzy vote to gene from patient)

	uniq_subnets = []
	count_hbeta = 0
	f_matrix = []
	print len(subnets)
	
	for s_name in subnets:
	
		S_nodes = subnets[s_name]
		
		# accounting for heterogeneity in both control and test
		if 'PFS' in option:
			beta_1 = {g: np.mean(weights_1[g]) for g in S_nodes}
			beta_2 = {g: np.mean(weights_2[g]) for g in S_nodes}
		else: 
		# SPS
			if symm == 'Y':
				snet_scores_1 = np.array([sum(weights_1[g][i] for g in S_nodes) for i in range(num_samples_1)])
				chosen_1 = list(np.argsort(-np.array(snet_scores_1))[0:k])
				beta_1 = {g: np.mean([weights_1[g][j] for j in chosen_1]) for g in S_nodes}
			else:
				beta_1 = {g: np.mean(weights_1[g]) for g in S_nodes}	
			
			if (len(weights_2[S_nodes[0]]) > k):
				snet_scores_2 = np.array([sum(weights_2[g][i] for g in S_nodes) for i in range(num_samples_2)])
				chosen_2 = list(np.argsort(-np.array(snet_scores_2))[0:k])
				beta_2 = {g: np.mean([weights_2[g][j] for j in chosen_2]) for g in S_nodes}
			else:
				beta_2 = {g: np.mean(weights_2[g]) for g in S_nodes}
	
		for j in range(num_samples_1):
			scores_1[s_name].append(sum(beta_1[g] * weights_1[g][j] for g in S_nodes))
			scores_2[s_name].append(sum(beta_2[g] * weights_1[g][j] for g in S_nodes))
		for i in range(num_samples_2):
			scores_1[s_name].append(sum(beta_1[g] * weights_2[g][i] for g in S_nodes))
			scores_2[s_name].append(sum(beta_2[g] * weights_2[g][i] for g in S_nodes))
	
		g_string = '\t'.join(sorted(S_nodes))
		
		# filter out duplicate subnets
		if not(g_string in uniq_subnets):
		
			uniq_subnets.append(g_string)
	
			beta_t = {g: np.mean(weights_2[g]) for g in S_nodes}
			if len(filter(lambda g: beta_t[g] > beta, S_nodes)) >= 5:
			
				count_hbeta += 1

				s1 = np.array(scores_1[s_name])
				s2 = np.array(scores_2[s_name])
		
				if np.var(s1) and np.var(s2):
					pval = R.r['t.test'](R.FloatVector(s2),R.FloatVector(s1),paired=True,alternative='greater')[2][0]
					if pval <= P_VAL_CUT_OFF:
						subnet2pval[s_name] = pval
						f_matrix.append(s1)
						f_matrix.append(s2)
						sgnf_subnets.append(s_name + '__S1')
						sgnf_subnets.append(s_name + '__S2')
						
	print 'subnets satisfying beta condition: ' + str(count_hbeta)
	
	plot_PCA('results', f_matrix, option, classes_c, classes_t, sgnf_subnets)
				
	return subnet2pval
	
# this function generates subnetworks and calculates their scores over the actual dataset
def snet(Xc, Xt, X_genes, subnets, theta_1, theta_2, beta, k, option, classes_c, classes_t, symm):

	num_c, num_t = np.shape(Xc)[1], np.shape(Xt)[1]
	
	"""computing weights"""

	print 'generating subnetworks...'

	# get gene weights in 2 patient groups
	weights_c = compute_weights(Xc, X_genes, theta_1, theta_2)
	weights_t = compute_weights(Xt, X_genes, theta_1, theta_2)

	"""subnetwork scoring"""

	print 'scoring subnetworks...'

	# for each sample, get subnet scores
	subnet_scores = get_subnet_scores (subnets, weights_c, weights_t, num_c, num_t, beta, k, option, classes_c, classes_t, symm)

	return subnet_scores

# filters GE matrix to contain expression of common genes in identical order
def rearrange_gdata(Xc, Xt, Xc_genes, Xt_genes):

	Xc_genes, Xt_genes = list(Xc_genes), list(Xt_genes)
	X_genes = list(set(Xc_genes) & set(Xt_genes))
	genes_ci, genes_ti = [], []

	# get indices for genes common to both groups
	for gene in X_genes:
		genes_ci.append(Xc_genes.index(gene))
		genes_ti.append(Xt_genes.index(gene))
	
	# filter expression matrix to include common genes only
	Xc = Xc[genes_ci]
	Xt = Xt[genes_ti]

	return Xc, Xt, X_genes

# write results to file
def write_to_file (f_path, subnets, sgnf_subnets):

	d_path = f_path.split('/')[0]

	# get path of the directory where script is running
	curr_path = os.path.dirname(os.path.abspath(__file__)) + '/'

	# normalize path as per local OS
	d_path = curr_path + os.path.normpath(d_path)
	f_path = curr_path + os.path.normpath(f_path)

	# create directory if it does not exist
	if not os.path.isdir(d_path):
		os.makedirs(d_path)

	subnet_genes = set()
	fo = open(f_path, 'wb')
	fo.write('subnet\tp_value\tgenes\n')
	for snet in sgnf_subnets:
		s_genes = '\t'.join(subnets[snet])
		subnet_genes |= set(subnets[snet])
		fo.write(snet + '\t' + str(sgnf_subnets[snet]) + '\t' + s_genes + '\n')
	fo.close()
	
# get the sample labels from the info file
def get_classes(fsi_name, snames):

	s2l = {}
	classes = []
	fsi = open(fsi_name, 'rb')
	#header = fsi.readline()
	for line in fsi:
		row = line.strip().split('\t')
		s2l[row[0]] = row[1]
	fsi.close()
	
	for s in snames:
		if s in s2l:
			classes.append(s2l[s])
		else:
			print 'No info found for sample ' + s + ' in the sample info file ' + fsi_name + '.'
			exit()
	
	return classes

# main function
def main(argv):

	start_time = time()
	
	'''argument handling'''

	fe_c_name, fe_t_name, fs_i_name, fpw_name, theta_1, theta_2, beta, k, option, symm = handle_args(argv)
	
	'''loading files'''

	print 'loading input files...'

	# control/normal group
	Xc, Xc_genes, snames_c, skipped_c = load_file(fe_c_name)
	if (skipped_c):
		print 'Control group: missing data -- skipped %d lines' %(skipped_c)

	# test/disease group
	Xt, Xt_genes, snames_t, skipped_t = load_file(fe_t_name)
	if (skipped_t):
		print 'Test group: missing data -- skipped %d lines' %(skipped_t)
		
	# classes -- control and test
	classes_c = get_classes(fs_i_name, snames_c) #len(snames_c)*['Solid Tissue Normal']
	classes_t = get_classes(fs_i_name, snames_t)

	# removing genes not common to test and control
	Xc, Xt, X_genes = rearrange_gdata(Xc, Xt, Xc_genes, Xt_genes)
	
	# pathways
	pw2G = load_pathways(fpw_name, X_genes)
	print len(pw2G)
	
	# form complexes
	print 'forming subnetworks...'
	subnets = break_into_subnets(pw2G)
	print len(subnets)
	
	# pfsnet
	sgnf_sn = snet(Xc, Xt, X_genes, subnets, theta_1, theta_2, beta, k, option, classes_c, classes_t, symm)

	# print results
	
	root = 'results/'
	write_to_file(root + option.upper() + '/' + option.lower() + '_b_' + str(beta) + '.txt', subnets, sgnf_sn)
	
	print 'DONE. Results can be found in the ' + root + ' directory.'

	print 'time elapsed (in seconds):', (time() - start_time)

	return

if __name__ == "__main__":
	main(sys.argv[1:])

