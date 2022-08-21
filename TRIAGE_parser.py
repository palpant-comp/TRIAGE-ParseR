###################################
### TRIAGE-ParseR               ###
###   Woo Jun (Chris) Shim      ###
###   contact: w.shim@uq.edu.au ###  
###################################

### Modules
import numpy as np
import pandas as pd
import os, sys, scipy.stats, collections, requests
import requests, json
from optparse import OptionParser
from sklearn.mixture import GaussianMixture

### Functions
def read_file_items(filename, col=0, numeric=False):
	results = []
	temp = open(filename, 'r')
	for i in temp:
		try:
			if not i[0].startswith('#'):
				i = i.strip().split()
				if numeric==False:
					results.append(i[col])
				else:
					results.append(float(i[col]))
		except:
			continue
	return results

def intersection(data1, data2):
    # takes a pair of lists and returns a list of shared items
    # data1 or data2 = a list of elements
    a = collections.Counter(data1)
    b = collections.Counter(data2)
    return list((a & b).elements())

def most_frequent_element(input_list): 
	# find an element in the list with the highest count
    counter = 0
    num = input_list[0]      
    for i in set(input_list): 
        curr_frequency = input_list.count(i) 
        if curr_frequency > counter: 
            counter = curr_frequency 
            num = i   
    return num

def write_table(input_data, output_file, numeric=False):
    # exports a table text file 
    # takes a dictionary of dictionaries (e.g. input_data[row1][col1]=xx, input_data[row1][col2]=xx, ..)    
    output_ = open(output_file, 'w')
    rownames = []
    for row in input_data:
        colnames = input_data[row].keys()
        rownames.append(row)
    if numeric==False:
        colnames = [col for col in colnames]
    else:
        colnames = [int(col) for col in colnames]
    colnames.sort()
    rownames.sort()
    first_line = ''
    for col in colnames:
        first_line += '\t'+str(col)
    output_.write(first_line+'\n')
    for row in rownames:
        line = str(row)
        for col in colnames:
            line += '\t'+str(input_data[row][str(col)])
        output_.write(line+'\n')

def find_gmm_cluster(input_list):
	temp = pd.read_csv(input_list, index_col='GENE')
	group = {}	
	for i in temp.columns:
		group[i] = []
	for j in range(temp.shape[0]):	# find the cluster for all genes
		id_ = temp.columns[np.argmax(temp.iloc[j])]
		group[id_].append(temp.index[j])
	return group

def create_sig_go_gmm_table(groups, ppi_threshold=0.01, enrichment_threshold=0.01, species=9606):
	# test ppi enrichment and if significant create a summary data (description, fdr, cluster)
	output_ = {}
	clusters = []
	for i in range(1, len(groups)+1):
		test = api_string(genes=groups['cluster'+str(i)], method='ppi_enrichment', output_format="tsv", species=species)
		#print (float(test[-1][-1]))
		if not test[0][0] =='Error':
			if float(test[-1][-1]) < ppi_threshold:  
				clusters.append(i)	
				if options.verbose==1:
					print ('cluster'+str(i),', PPI enrichment p-value =',float(test[-1][-1]))
	#print (ppi_threshold, clusters)
	for i in clusters:
		results = api_string(genes=groups['cluster'+str(i)], method='enrichment', species=species, enrichment_threshold=enrichment_threshold)	
		output__ = [[results[i][-1], float(results[i][-2])] for i in range(len(results))]		
		for j in output__:
			name = j[0].replace(' ','_')		
			if name not in output_:
				output_[name] = {}
				for m in clusters:
					output_[name]['cluster'+str(m)] = 1
			output_[name]['cluster'+str(i)] = j[1]
	return output_

def api_string(genes, method='ppi_enrichment', output_format="json", species=9606, enrichment_threshold=0.01):
	# returns enrichment data from STRING using API
	# GO enrichment, method='enrichment'
	# PPI enrichment, method='ppi_enrichment'
	string_api_url = "https://version-11-5.string-db.org/api"
	request_url = "/".join([string_api_url, output_format, method])
	params = {"identifiers" : "%0d".join(genes), "species" : species, "caller_identity" : "comodulation"}
	response = requests.post(request_url, data=params)
	results = []
	if method=='ppi_enrichment':	# tsv 
		for line in response.text.strip().split("\n"):
			results.append(line.split("\t"))
	else:
		data = json.loads(response.text)	# json
		for row in data:
			term = row["term"]
			preferred_names = ",".join(row["preferredNames"])
			fdr = float(row["fdr"])
			description = row["description"]
			category = row["category"]
			if category == "Process" and fdr < 0.01:
				results.append([term, preferred_names, str(fdr), description])    
	return results

def write_file(data, filename):	
	temp = open(filename, 'w')
	line = '#gene'
	for item in data.columns:
		line += '\t'+str(item)
	temp.write(line+'\n')
	genes = data.index
	for i in range(data.shape[0]):
		line = genes[i]
		for j in range(data.shape[1]):
			line += '\t'+str(data.iloc[i,j])
		temp.write(line+'\n')
	temp.close

def calculate_variance(input):
	return [np.var(input.iloc[:,i]) for i in range(input.shape[1])]


if __name__ == '__main__':
    ### Command line options   
	parser = OptionParser()
	parser.add_option('-i', '--i', dest='input', help='input')    
	parser.add_option('-r', '--r', dest='H3K27me3_pc', help='pre-calculated H3K27me3 principal components', default='./data/pca_roadmap')
	parser.add_option('-p', '--p', dest='no_pca', help='Number of PCs to use', default=10, type=int)
	parser.add_option('-g', '--g', dest='no_genes', help='Number of top genes to use', default=100, type=int)
	parser.add_option('-t', '--t', dest='no_iter', help='Number of iterations for model selection', default=100, type=int)
	parser.add_option('-o', '--o', dest='output_directory', help='Output directory', default='./results/')
	parser.add_option('-e', '--e', dest='go_analysis', help='Whether to perform GO enrichment analysis (1: Yes, default, 0: No', default=1, type=int)
	parser.add_option('-a', '--a', dest='input_type', help='Input type (option: table or list)', default='list')
	parser.add_option('-v', '--v', dest='verbose', help='Level of verbose (option: 1 or 0', default=1)
	parser.add_option('-w', '--w', dest='max_cluster', help='Max. number of clusters', default=10)
	parser.add_option('-j', '--j', dest='gene_order', help='Gene sort direction (option: ascending or descending)', default='descending')

	options = parser.parse_args()[0]  

	if options.input == None:
		sys.exit('Exiting: input (TRIAGE output table or a list of genes) is required (-i)')
	else:
		if not os.path.exists(options.output_directory):
			os.mkdir(options.output_directory)
		if not os.path.exists(options.output_directory+'/gene_clusters'):
			os.mkdir(options.output_directory+'/gene_clusters')

		### 1. Load data
		pca = pd.read_csv(options.H3K27me3_pc+'.csv', index_col='GENE')
		rows = pca.index
		cols = pca.head()
		if options.input_type=='table':
			disc = pd.read_csv(options.input, index_col='GENE', sep='\t')
			labels = disc.head()
		else:
			labels = ['output']
		if options.input_type=='list':
			genes = read_file_items(options.input)
			options.no_genes = len(genes)
		if options.verbose==1:			
			print ("Using", options.no_pca, 'top PCs from top', options.no_genes, 'genes with',options.no_iter,'iterations for model selection')

		### 2. Find gene clusters
		if options.verbose==1:	
			print ('Finding gene clusters ...')
		for label in labels:
			if options.verbose==1:
				print (label)
			if options.input_type=='table':
				if options.gene_order=='descending':	
					genes = disc.sort_values(by=[label], ascending=False).index[:options.no_genes]
				else:
					genes = disc.sort_values(by=[label], ascending=True).index[:options.no_genes]
			# else:
			# 	genes = read_file_items(options.input)
			genes = intersection(genes, rows)
			aa = pca.loc[genes]

			# find PCs with most variance among the gene set
			cc = calculate_variance(aa)
			idx = list(np.argsort(cc))
			idx.reverse()
			aa = aa.iloc[:,idx[:options.no_pca]]

			# find an optimal number of clusters using Bayesian information criterion 
			bic_ = []
			for no in range(options.no_iter):
				bic_result = []
				for n in range(2, options.max_cluster):
					gmm = GaussianMixture(n).fit(aa)  
					bic_result.append(gmm.bic(aa))
				bic_.append(bic_result.index(np.min(bic_result)) + 1)  
			n = most_frequent_element(input_list=bic_)
			if options.verbose==1:
				print ('Number of cluster =', n)
			gmm = GaussianMixture(n, random_state=42).fit(aa)
			bb = pd.DataFrame(gmm.predict_proba(aa))
			bb.index = genes
			results = pd.DataFrame(columns=['cluster'+str(i+1) for i in range(n)], index=genes)

			for i in range(len(genes)):
				results.loc[genes[i]] = list(bb.iloc[i])
			results.to_csv(options.output_directory+'/gene_clusters/'+label+'_gene_clusters.txt', index_label='GENE')

		### 3. GO ENRICHMENT 
		if options.go_analysis==1:
			if options.verbose==1:
				print ('Performing GO enrichment analysis ... ')
			if not os.path.exists(options.output_directory+'/go'):
				os.mkdir(options.output_directory+'/go')
			files = os.listdir(options.output_directory+'/gene_clusters/')			
			for i in files:
				a = i.split('_gene_clusters')
				name = a[0]
				if options.verbose==1:
					print (a[0])
				groups = find_gmm_cluster(input_list=options.output_directory+'/gene_clusters/'+name+'_gene_clusters.txt')
				#print (groups)
				results = create_sig_go_gmm_table(groups=groups, ppi_threshold=0.01, enrichment_threshold=0.01)
				write_table(results, options.output_directory+'/go/'+name+'_go.txt')
