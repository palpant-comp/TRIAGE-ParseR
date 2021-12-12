### PCA-GMM FOR TRIAGE GENE CLUSTERING ###
### DEVELOPED AND WRITTEN BY WOO JUN SHIM ###

import numpy as np
import os, sys, scipy.stats, collections, requests
from scipy.stats import beta
from optparse import OptionParser
from sklearn.mixture import GaussianMixture

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

def read_table(input_file, numeric=True, as_list=False, include_sum=False, threshold_=None, del_colnames=True):
    ### reads in a table file and returns a dictionary
    ### Note. the number of colnames must be always n-1 where n is the total number of columns (including rownames)
    ### numerical = a boolean to indicate whether the data should be read as numbers (True) or string (False)
    ### as_list = a boolean to indicate whether the data is to be stored in the form of dictionary of lists, rather than dictionary of dictionaries
    ### include_sum = a boolean to indicate whether total counts for each row is to be recorded (only works for numerical data)
    ### threshold = a numerical value to set a threshold value entries with values below which are not included (if threshold_= None --> no threshold)
    ### del_colnames = a boolean to indicate whether to remove 'colnames' in the table
    ### e.g. output[row][col] = xx (when as_list=False) OR output[row] = [xx, xx, ..] (when as_list=True)
    ###      output[row]['total'] = sum of total for the row
    ###      output['colnames'] = [a list of column names]
    results = {}
    results['colnames'] = []
    input_ = open(input_file, 'r')
    cnt= 0
    for line in input_:
        line = line.replace(' ','_')
        line = line.strip().split()
        sum_ = 0.0
        if cnt == 0:
            for item in line:
                results['colnames'].append(item)
            cnt = 1
        else:
            if line[0] not in results:
                if as_list == False:
                    results[line[0]] = {}
                else:
                    results[line[0]] = []
            for no in range(len(results['colnames'])):
                group_name = results['colnames'][no]                
                if numeric == True:
                    if as_list == False:
                        if threshold_!= None:
                            if float(line[no+1]) >= float(threshold_):
                                results[line[0]][group_name] = float(line[no+1])
                        else:
                            results[line[0]][group_name] = float(line[no + 1])
                    else:
                        if threshold_!= None:
                            if float(line[no + 1]) >= float(threshold_):
                                results[line[0]].append(float(line[no + 1]))
                        else:
                            results[line[0]].append([group_name, float(line[no+1])])
                    if include_sum == True:
                        sum_ += float(line[no+1])
                else:
                    if as_list == False:
                        results[line[0]][group_name] = line[no+1]
                    else:
                        results[line[0]].append(line[no+1])
            if include_sum == True:
                results[line[0]]['sum'] = sum_
    if del_colnames==True:
        del results['colnames']
    return results

def find_index_for_array(a,b):
	### a = a subset, b = a superset set
	a = np.array(a)
	b = np.array(b)
	sorter = np.argsort(b)
	return sorter[np.searchsorted(b, a, sorter=sorter)]

def initiate_table(rows, cols, default=0):
    results = {}
    for r in rows:
        results[r] = {}
        for c in cols:
            results[r][c] = default
    return results

def get_top_rownames(input_table, col, top_no=100, reverse=True):
	### returns rownames of a table with the top values in a column
	input_table = check_table(input_table)
	temp = [[g, float(input_table[g][col])] for g in input_table]
	temp = order_list(temp, col=1, reverse=reverse)
	if top_no==None:
		top_no = len(get_rownames(temp))
	return [i[0] for i in temp[:top_no]]

def order_list(input_list, col, reverse = False):
    return sorted(input_list, key=lambda x: float(x[col]), reverse=reverse)

def check_table(filename, numeric=True):
	if type(filename)==str:
		filename = read_table(filename, numeric=numeric)
	return filename

def intersection(data1, data2):
    ''' takes a pair of lists and returns a set of intersected items
    data1 or data2 = a list of elements
    '''
    a = collections.Counter(data1)
    b = collections.Counter(data2)
    return list((a & b).elements())

def most_frequent_element(input_list): 
    counter = 0
    num = input_list[0] 
      
    for i in input_list: 
        curr_frequency = input_list.count(i) 
        if curr_frequency > counter: 
            counter = curr_frequency 
            num = i   
    return num

def write_table(input_data, output_file, numeric=False):
    ''' writes out a table text file 
    takes a dictionary of dictionaries (e.g. input_data[row1][col1]=xx, input_data[row1][col2]=xx, ..) 
    '''
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
	### input_list is an output from PCA-GMM cluster
	### returns a dictionary of lists
	temp = check_list(input_list)
	group = {}
	for i in range(1,len(temp[0])):
		group[i] = []
	for j in range(len(temp)):	# find the cluster for all genes 
		a = [float(temp[j][m]) for m in range(1, len(temp[j]))]
		idx = a.index(np.max(a))+1
		group[idx].append(temp[j][0])
	return group

def create_sig_go_gmm_table(groups, ppi_threshold=0.001, top_no=10):
	### test ppi enrichment and if significant create a summary data (description, fdr, cluster)
	output_ = {}
	clusters = []
	for i in range(1, len(groups)+1):
		test = api_string(genes=groups[i], method='ppi_enrichment', output_format="tsv")
		#print (float(test[-1][-1]))
		if float(test[-1][-1]) < ppi_threshold:  
			clusters.append(i)	
	for i in clusters:
		results = api_string(genes=groups[i], method='enrichment')			
		for j in results[:top_no]:
			name = j[-1].replace(' ','_')		
			if name not in output_:
				output_[name] = {}
				for m in clusters:
					output_[name]['cluster'+str(m)] = 1
			output_[name]['cluster'+str(i)] = j[-2]
	return output_

def api_string(genes, method='ppi_enrichment', output_format="tsv-no-header"):
	### returns enrichment data from STRING using API
	### GO enrichment, method='enrichment'
	### PPI enrichment, method='ppi_enrichment'
	string_api_url = "https://string-db.org/api"
	request_url = "/".join([string_api_url, output_format, method])
	params = {"identifiers" : "%0d".join(genes), "species" : 9606, "caller_identity" : "pca_gmm"}
	response = requests.post(request_url, data=params)
	results = []
	for line in response.text.strip().split("\n"):
		results.append(line.split("\t"))
	return results

def get_colnames(table):
    results = []
    for gene in table:
        if len(results)==0:
            for item in table[gene]:
                results.append(item)
        else:
            break
    return results

def write_file(data, filename):
	temp = open(filename, 'w')
	if type(data)==dict:
		for key in data:
			line = str(key)
			if len(data[key]) != 0:
				for j in data[key]:
					line += '\t'+str(j)
			else:
				line += '\t'+str(data[key])
			line += '\n'
			temp.write(line)
		temp.close()
	else:
		for i in data:
			if len(i) != 0:
				line = str(i[0])
				for j in i[1:]:
					line += '\t'+str(j)
				line += '\n'
			temp.write(line)
		temp.close()

def check_list(input_file):
	if type(input_file)==str:
		results = read_file_as_list(input_file)
	else:
		results = input_file
	return results

def read_file_as_list(filename, col=None, numeric_col=None, keep_header=False):
	results = []
	temp = open(filename, 'r')
	for i in temp:
		i = i.replace(' ','_')
		i = i.strip().split()
		if len(i) > 0:
			if i[0].startswith('#'):
				if keep_header==True:
					header = i
				else:
					header = []	
				continue
			else:
				results.append([])
				if col==None:
					for no in range(len(i)):
						c = i[no]
						if numeric_col!=None:
							if no in numeric_col:						
								c = float(i[no])
						results[-1].extend([c])
				else:
					for c in col:
						if numeric_col!=None:
							if c in numeric_col:
								value = float(i[c])
						else:
							value = i[c]
						results[-1].extend(value)
	if keep_header==False:
		return results
	else:
		return (results, header)

if __name__ == '__main__':

    ### Command line options   
	parser = OptionParser()
	parser.add_option('-i', '--i', dest='discordance_table', help='input discordance table')    
	parser.add_option('-r', '--r', dest='H3K27me3_pc', help='pre-calculated H3K27me3 principal components', default='./data/pca_x')
	parser.add_option('-p', '--p', dest='no_pca', help='Number of PCs to use', default=67, type=int)
	parser.add_option('-g', '--g', dest='no_genes', help='Number of top genes to use', default=100, type=int)
	parser.add_option('-t', '--t', dest='no_iter', help='Number of iterations for model selection', default=10, type=int)
	parser.add_option('-o', '--o', dest='output_directory', help='Output directory', default='./results/')
	parser.add_option('-e', '--e', dest='go_analysis', help='Whether to perform GO enrichment analysis (1: Yes, default, 0: No', default=1, type=int)

	options = parser.parse_args()[0]  

	if options.discordance_table == None:
		sys.exit('Exiting: input discordance table is required (-i)')
	else:
		if not os.path.exists(options.output_directory):
			os.mkdir(options.output_directory)
		if not os.path.exists(options.output_directory+'/gene_clusters'):
			os.mkdir(options.output_directory+'/gene_clusters')
		no_pca = options.no_pca
		no_genes = options.no_genes
		no_iter = options.no_iter

		### 1. Load data
		pca = np.loadtxt(options.H3K27me3_pc+'.csv',delimiter=',', dtype=float)
		cols = read_file_items(options.H3K27me3_pc+'_cols.txt')
		rows = read_file_items(options.H3K27me3_pc+'_rows.txt')
		disc = read_table(options.discordance_table)
		col_idx = find_index_for_array(a=['PC'+str(i+1) for i in range(no_pca)],b=cols)
		labels = get_colnames(disc)
		results = initiate_table(rows=['PC'+str(i+1) for i in range(no_pca)], cols=labels, default=1)
		total_rows = len(rows)
		print ("H3K27me3 PCA table =",pca.shape)
		print ("Using", no_pca, 'PCs from top', no_genes, 'genes with',no_iter,'iterations for model selection')

		### 2. Find gene clusters
		# n_comp = [i+1 for i in range(no_pca)]		
		# print ('Finding gene clusters ...')
		# for label in labels:			
		# 	genes = get_top_rownames(input_table=disc, col=label, top_no=2000, reverse=True)
		# 	genes = intersection(genes, rows) 
		# 	row_idx = find_index_for_array(a=genes,b=rows)[:no_genes]  
		# 	aa = pca[row_idx,:]
		# 	bic_ = []
		# 	for no in range(no_iter):
		# 		bic_result = []
		# 		for n in n_comp:
		# 			gmm = GaussianMixture(n).fit(aa)  
		# 			bic_result.append(gmm.bic(aa))
		# 		bic_.append(bic_result.index(np.min(bic_result)) + 1)  # optimal number of pcincipal components 
		# 	n = most_frequent_element(input_list=bic_)
		# 	print (label, n)
		# 	gmm = GaussianMixture(n, random_state=42).fit(aa)
		# 	bb = gmm.predict_proba(aa)  # Assess top x TRIAGE-prioritised genes
		# 	results = [['#gene']+['cluster'+str(i+1) for i in range(n)]]
		# 	for i in range(len(bb)):
		# 		gene = rows[row_idx[i]]
		# 		results.append([gene])
		# 		for j in range(len(bb[i])):
		# 			results[-1].extend([bb[i][j]])
		# 	write_file(results, options.output_directory+'/gene_clusters/'+label+'_gene_clusters.txt')

		# 3. GO ENRICHMENT 
		if options.go_analysis==1:
			print ('Performing GO enrichment analysis ... ')
			if not os.path.exists(options.output_directory+'/go'):
				os.mkdir(options.output_directory+'/go')
			files = os.listdir(options.output_directory+'/gene_clusters/')			
			for i in files:
				try:
					a = i.split('_gene_clusters')
					name = a[0]
					print (a[0])
					groups = find_gmm_cluster(input_list=options.output_directory+'/gene_clusters/'+name+'_gene_clusters.txt')
					results = create_sig_go_gmm_table(groups=groups, ppi_threshold=0.001, top_no=10)			
					write_table(results, options.output_directory+'/go/'+name+'_go.txt')
				except:
					continue

