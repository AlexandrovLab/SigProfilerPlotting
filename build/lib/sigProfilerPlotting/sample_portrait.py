#!/usr/bin/env python3
 
#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu


import re
import os
import sys
import argparse
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mplpatches



def samplePortrait (sample_matrices_path, output_path, project, percentage=False):
	pcawg = False
	sig_probs = False
	if sample_matrices_path[-1] != "/":
		sample_matrices_path += "/"
	sample_matrices_path += "output/"


	pp = PdfPages(output_path + 'sample_portrait_' + project + '.pdf')

	mutations_6 = OrderedDict()
	mutations_24 = OrderedDict()
	mutations_96 = OrderedDict()
	mutations_384 = OrderedDict()
	mutations_1536 = OrderedDict()
	mutations_6144 = OrderedDict()
	mutations_5 = OrderedDict()
	mutations_3 = OrderedDict()

	mutations_78 = OrderedDict()
	mutations_312 = OrderedDict()

	mutations_83 = OrderedDict()
	mutations_ID96 = OrderedDict()
	mutations_simple = OrderedDict()


	file_example = os.listdir(sample_matrices_path + "SBS/")[1]
	# if '.DS_Store' in file_example:
	# 	file_example.remove('.DS_Store')
	file_extension = file_example.split(".")[-1]
	file_extension = "all"
	SBS96 = True
	SBS6 = True
	SBS24 = True
	SBS384 = True
	SBS1536 = True
	SBS6144 = True

	ID28 = True
	ID83 = True
	ID96 = True

	DBS78 = True
	DBS312 = True


	########### Gather mutation counts ##############################
	try:
		samples = []
		with open(sample_matrices_path + "SBS/" + project + ".SBS96." + file_extension) as f:
			next(f)
			first_line = f.readline()
			first_line = first_line.strip().split()
			if first_line[0][1] == ">":
				pcawg = True
			if first_line[0][5] != "]" and first_line[0][1] != ">":
				sys.exit("The matrix does not match the correct SBS96 format. Please check you formatting and rerun this plotting function.")
			total_count = []

		with open (sample_matrices_path + "SBS/" + project + ".SBS96." + file_extension) as f:
			first_line = f.readline()
			if pcawg:
				samples = first_line.strip().split(",")
				samples = samples[2:]
			else:
				samples = first_line.strip().split("\t")
				samples = samples[1:]
				
			for sample in samples:
				mutations_96[sample] = OrderedDict()
				mutations_96[sample]['C>A'] = OrderedDict()
				mutations_96[sample]['C>G'] = OrderedDict()
				mutations_96[sample]['C>T'] = OrderedDict()
				mutations_96[sample]['T>A'] = OrderedDict()
				mutations_96[sample]['T>C'] = OrderedDict()
				mutations_96[sample]['T>G'] = OrderedDict()

			for lines in f:
				if pcawg:
					line = lines.strip().split(",")
					mut_type = line[0]
					nuc = line[1][0] + "[" + mut_type + "]" + line[1][2]
					sample_index = 2
				else:
					line = lines.strip().split()
					nuc = line[0]
					mut_type = line[0][2:5]
					sample_index = 1

				for sample in samples:
					if percentage:
						mutCount = float(line[sample_index])
						if mutCount < 1 and mutCount > 0:
							sig_probs = True
					else:
						mutCount = int(line[sample_index])
					mutations_96[sample][mut_type][nuc] = mutCount
					sample_index += 1
	except:
		SBS96 = False
		print("No SBS-96 provided")


	########################### SBS-6 ################################################
	if True:
	#try:
		with open (sample_matrices_path + "SBS/" + project + ".SBS6." + file_extension) as f:
			first_line = f.readline()
			for sample in samples:
				mutations_6[sample] = OrderedDict()
				mutations_6[sample]['C>A'] = 0
				mutations_6[sample]['C>G'] = 0
				mutations_6[sample]['C>T'] = 0
				mutations_6[sample]['T>A'] = 0
				mutations_6[sample]['T>C'] = 0
				mutations_6[sample]['T>G'] = 0

			for lines in f:
				line = lines.strip().split()
				nuc = line[0]
				mut_type = line[0]
				sample_index = 1

				for sample in samples:
					if percentage:
						mutCount = float(line[sample_index])
						if mutCount < 1 and mutCount > 0:
							sig_probs = True
					else:
						mutCount = int(line[sample_index])
					mutations_6[sample][mut_type] = mutCount
					sample_index += 1

	#except:
	#	SBS6 = False
	#	print("No SBS-6 provided")

		########################### SBS-24 ################################################
	try:
		with open (sample_matrices_path + "SBS/" + project + ".SBS24." + file_extension) as f:
			first_line = f.readline()
			for sample in samples:
				mutations_24[sample] = OrderedDict()
				mutations_24[sample]['C>A'] = [0,0]
				mutations_24[sample]['C>G'] = [0,0]
				mutations_24[sample]['C>T'] = [0,0]
				mutations_24[sample]['T>A'] = [0,0]
				mutations_24[sample]['T>C'] = [0,0]
				mutations_24[sample]['T>G'] = [0,0]

			for lines in f:
				line = lines.strip().split()
				nuc = line[0][2:]
				bias = line[0][0]
				if bias == 'N' or bias == 'B':
					continue
				else:
					sample_index = 1
					for sample in samples:
						if percentage:
							mutCount = float(line[sample_index])
							if mutCount < 1 and mutCount > 0:
								sig_probs = True
						else:
							mutCount = int(line[sample_index])						
						if bias == 'T':
							mutations_24[sample][nuc][0] = mutCount
						else:
							mutations_24[sample][nuc][1] = mutCount
						sample_index += 1

	except:
		SBS24 = False
		print("No SBS-24 provided")



	########################### SBS-384 ################################################
	try:
		with open (sample_matrices_path + "SBS/" + project + ".SBS384." + file_extension) as f:
			first_line = f.readline()
			for sample in samples:
				mutations_384[sample] = OrderedDict()
				mutations_384[sample]['C>A'] = OrderedDict()
				mutations_384[sample]['C>G'] = OrderedDict()
				mutations_384[sample]['C>T'] = OrderedDict()
				mutations_384[sample]['T>A'] = OrderedDict()
				mutations_384[sample]['T>C'] = OrderedDict()
				mutations_384[sample]['T>G'] = OrderedDict()


			for lines in f:
				if pcawg:
					line = lines.strip().split(",")
					line = [x.replace('"','') for x in line]
					nuc = line[2][0] + "[" + line[1] + "]" + line[2][2]
					bias = line[0][0]
				else:
					line = lines.strip().split()
					nuc = line[0][2:]
					bias = line[0][0]
				if bias == 'N' or bias == 'B':
					continue
				else:
					if pcawg:
						mut_type = line[1]
						sample_index = 3
					else:					
						mut_type = line[0][4:7]
						sample_index = 1

					for sample in samples:
						if percentage:
							mutCount = float(line[sample_index])
							if mutCount < 1 and mutCount > 0:
								sig_probs = True
						else:
							mutCount = int(line[sample_index])
						if nuc not in mutations_384[sample][mut_type].keys():
							mutations_384[sample][mut_type][nuc] = [0,0]
						if bias == 'T':
							mutations_384[sample][mut_type][nuc][0] = mutCount
						else:
							mutations_384[sample][mut_type][nuc][1] = mutCount
						sample_index += 1

	except:
		SBS384 = False
		print("No SBS-384 provided")		

	########################### SBS-1536 ################################################
	try:
		mutations_1536 = OrderedDict()
		mutations_5 = OrderedDict()
		mutations_3 = OrderedDict()
		max_count = {}
		max_all = {}
		max_5 = {}
		max_3 = {}
		total_count = []
		total_counts = {'TT':0, 'TG':0,'TC':0,'TA':0,
						'GT':0,'GG':0,'GC':0,'GA':0,
						'CT':0,'CG':0,'CC':0,'CA':0,
						'AT':0,'AG':0,'AC':0,'AA':0,}
		total_counts_5 = {'T':0, 'G':0,'C':0,'A':0}
		total_counts_3 = {'T':0, 'G':0,'C':0,'A':0}

		with open (sample_matrices_path + "SBS/" + project + ".SBS1536." + file_extension) as f:
			first_line = f.readline()
			if pcawg:
				samples = first_line.strip().split(",")
				samples = samples[2:]
				samples = [x.replace('"','') for x in samples]
			# else:
			# 	samples = first_line.strip().split("\t")
			# 	samples = samples[1:]

			for sample in samples:
				max_all[sample] = 0
				max_5[sample] = 0
				max_3[sample] = 0
				total_counts[sample]= {'TT':0, 'TG':0,'TC':0,'TA':0,
									   'GT':0,'GG':0,'GC':0,'GA':0,
									   'CT':0,'CG':0,'CC':0,'CA':0,
									   'AT':0,'AG':0,'AC':0,'AA':0,}
				total_counts_5[sample]= {'T':0, 'G':0,'C':0,'A':0}
				total_counts_3[sample]= {'T':0, 'G':0,'C':0,'A':0}

				# mutations_96[sample] = OrderedDict()
				# mutations_96[sample]['C>A'] = OrderedDict()
				# mutations_96[sample]['C>G'] = OrderedDict()
				# mutations_96[sample]['C>T'] = OrderedDict()
				# mutations_96[sample]['T>A'] = OrderedDict()
				# mutations_96[sample]['T>C'] = OrderedDict()
				# mutations_96[sample]['T>G'] = OrderedDict()

				max_count[sample] = 0
				mutations_1536[sample] = OrderedDict()
				mutations_5[sample] = OrderedDict()
				mutations_3[sample] = OrderedDict()

				mutations_1536[sample]['C>A'] = OrderedDict()
				mutations_5[sample]['C>A'] = OrderedDict()
				mutations_3[sample]['C>A'] = OrderedDict()
				mutations_1536[sample]['C>A'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
											'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
											'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(), 
											'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
				mutations_5[sample]['C>A'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
				mutations_3[sample]['C>A'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}

				mutations_1536[sample]['C>G'] = OrderedDict()
				mutations_5[sample]['C>G'] = OrderedDict()
				mutations_3[sample]['C>G'] = OrderedDict()
				mutations_1536[sample]['C>G'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
											'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
											'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(), 
											'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
				mutations_5[sample]['C>G'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
				mutations_3[sample]['C>G'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}

				mutations_1536[sample]['C>T'] = OrderedDict()
				mutations_5[sample]['C>T'] = OrderedDict()
				mutations_3[sample]['C>T'] = OrderedDict()
				mutations_1536[sample]['C>T'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
											'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
											'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(), 
											'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
				mutations_5[sample]['C>T'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
				mutations_3[sample]['C>T'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}

				mutations_1536[sample]['T>A'] = OrderedDict()
				mutations_5[sample]['T>A'] = OrderedDict()
				mutations_3[sample]['T>A'] = OrderedDict()
				mutations_1536[sample]['T>A'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
											'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
											'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(), 
											'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
				mutations_5[sample]['T>A'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
				mutations_3[sample]['T>A'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}

				mutations_1536[sample]['T>C'] = OrderedDict()
				mutations_5[sample]['T>C'] = OrderedDict()
				mutations_3[sample]['T>C'] = OrderedDict()
				mutations_1536[sample]['T>C'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
											'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
											'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(), 
											'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
				mutations_5[sample]['T>C'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
				mutations_3[sample]['T>C'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}

				mutations_1536[sample]['T>G'] = OrderedDict()
				mutations_5[sample]['T>G'] = OrderedDict()
				mutations_3[sample]['T>G'] = OrderedDict()
				mutations_1536[sample]['T>G'] = {'TT':OrderedDict(), 'TG':OrderedDict(), 'TC':OrderedDict(), 'TA':OrderedDict(),
											'GT':OrderedDict(), 'GG':OrderedDict(), 'GC':OrderedDict(), 'GA':OrderedDict(),
											'CT':OrderedDict(), 'CG':OrderedDict(), 'CC':OrderedDict(), 'CA':OrderedDict(), 
											'AT':OrderedDict(), 'AG':OrderedDict(), 'AC':OrderedDict(), 'AA':OrderedDict()}
				mutations_5[sample]['T>G'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}
				mutations_3[sample]['T>G'] = {'T':OrderedDict(),'G':OrderedDict(),'C':OrderedDict(),'A':OrderedDict()}


			for lines in f:
				if pcawg:
					line = lines.strip().split(",")
					line = [x.replace('"','') for x in line]
					nuc = line[1][0:2] + "[" + line[0] + "]" + line[1][3:]
					mut_type = line[0]
					penta_key = line[1][0] + line[1][-1]
					tri_key = line[1][1] + line[1][-2]
					sample_index = 2		
				else:			
					line = lines.strip().split()
					nuc = line[0]
					mut_type = line[0][3:6]
					penta_key = line[0][0] + line[0][-1]
					tri_key = line[0][1] + line[0][-2]
					sample_index = 1

					tri = line[0][1:8]

				for sample in samples:

					if tri not in mutations_96[sample][mut_type]:
						mutations_96[sample][mut_type][tri] = 0
					if percentage:
						mutCount = float(line[sample_index])
						if mutCount < 1 and mutCount > 0:
							sig_probs = True
					else:
						mutCount = int(line[sample_index])

					if pcawg:
						sample_ref = sample_index - 2
					else:
						sample_ref = sample_index - 1
					if mutCount > max_count[samples[sample_ref]]:
						max_count[samples[sample_ref]] = mutCount

					if mutCount > max_all[sample]:
						max_all[sample] = mutCount

					mutations_1536[sample][mut_type][penta_key][tri_key] = mutCount
					total_counts[sample][penta_key] += mutCount
					total_counts_5[sample][penta_key[0]] += mutCount
					total_counts_3[sample][penta_key[1]] += mutCount
					penta_key_short = penta_key[0]
					mutations_5[sample][mut_type][penta_key_short][tri_key] = 0
					mutations_3[sample][mut_type][penta_key_short][tri_key] = 0
					#mutations_96[sample][mut_type][tri] += mutCount
					sample_index += 1

	except:
		SBS1536 = False
		print("No SBS-1536 provided")		


	########################### ID-83 ################################################
	try:
		indel_types = ['1:Del:C:0', '1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5',
					   '1:Del:T:0', '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5',
					   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
					   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5', 
							# >1bp INDELS
					   '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5',
					   '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5',
					   '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5',
					   '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5',
					   '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5', 
					   '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5', 
					   '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5',
					   '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5',
							#MicroHomology INDELS
					   '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3',
					   '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5', '2:Ins:M:1', 
					   '3:Ins:M:1', '3:Ins:M:2', '4:Ins:M:1', '4:Ins:M:2', '4:Ins:M:3', '5:Ins:M:1', 
					   '5:Ins:M:2', '5:Ins:M:3', '5:Ins:M:4', '5:Ins:M:5']

		mutations = OrderedDict()
		with open (sample_matrices_path + "ID/" + project + ".ID83." + file_extension) as f:
			first_line = f.readline()
			for sample in samples:
				mutations_83[sample] = OrderedDict()
				mutations_83[sample]['1DelC'] = [0,0,0,0,0,0]
				mutations_83[sample]['1DelT'] = [0,0,0,0,0,0]
				mutations_83[sample]['1InsC'] = [0,0,0,0,0,0]
				mutations_83[sample]['1InsT'] = [0,0,0,0,0,0]
				mutations_83[sample]['2DelR'] = [0,0,0,0,0,0]
				mutations_83[sample]['3DelR'] = [0,0,0,0,0,0]
				mutations_83[sample]['4DelR'] = [0,0,0,0,0,0]
				mutations_83[sample]['5DelR'] = [0,0,0,0,0,0]
				mutations_83[sample]['2InsR'] = [0,0,0,0,0,0]
				mutations_83[sample]['3InsR'] = [0,0,0,0,0,0]
				mutations_83[sample]['3InsR'] = [0,0,0,0,0,0]
				mutations_83[sample]['4InsR'] = [0,0,0,0,0,0]
				mutations_83[sample]['5InsR'] = [0,0,0,0,0,0]
				mutations_83[sample]['2DelM'] = [0]
				mutations_83[sample]['3DelM'] = [0,0]
				mutations_83[sample]['4DelM'] = [0,0,0]
				mutations_83[sample]['5DelM'] = [0,0,0,0,0]

			for lines in f:
				if pcawg:
					line = lines.strip().split(",")
					line = [x.replace('"','') for x in line]
					if line[1] == 'repeats':
						mut_type = line[2][0] + line[0][0] + line[0][1].lower() + line[0][2].lower() + "R"
					else:
						mut_type = line[2][0] + line[0][0] + line[0][1].lower() + line[0][2].lower() + line[1][0]
					try:
						repeat_size = int(line[3])
					except:
						repeat_size = int(line[3][0])
					if line[1] == 'MH':
						repeat_size -= 1
					sample_index = 4
				else:
					line = lines.strip().split()
					if line[0] not in indel_types:
						continue
					categories = line[0].split(":")
					mut_type = categories[0] + categories[1] + categories[2]
					repeat_size = int(categories[3])
					if categories[2] == 'M':
						repeat_size -= 1
					sample_index = 1

				for sample in samples:
					if mut_type in mutations_83[sample].keys():
						if percentage:
							mutCount = float(line[sample_index])
							if mutCount < 1 and mutCount > 0:
								sig_probs = True
						else:
							mutCount = int(line[sample_index])
						mutations_83[sample][mut_type][repeat_size] = mutCount
					else:
						continue
					sample_index += 1

	except:
		ID83 = False
		print("No ID83 provided")



########################### ID-96 ################################################
	try:
		indel_types_tsb = []
		tsb_I = ['T','U','N','B','Q']
		indel_types = ['1:Del:C:0', '1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5',
					   '1:Del:T:0', '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5',
					   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
					   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5', 
							# >1bp INDELS
					   '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5',
					   '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5',
					   '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5',
					   '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5',
					   '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5', 
					   '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5', 
					   '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5',
					   '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5',
							#MicroHomology INDELS
					   '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3',
					   '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5']

		for indels in indel_types:
			for tsbs in tsb_I:
				indel_types_tsb.append(tsbs + ":" + indels)					   

		sig_probs = False
		with open (sample_matrices_path + "ID/" + project + ".ID415." + file_extension) as f:
			first_line = f.readline()
			if pcawg:
				samples = first_line.strip().split(",")
				samples = samples[4:]
				samples = [x.replace('"','') for x in samples]
			else:
				pass			
				#samples = first_line.strip().split("\t")
				#samples = samples[1:]
			for sample in samples:
				mutations_ID96[sample] = OrderedDict()
				mutations_ID96[sample]['1DelC'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['1DelT'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['1InsC'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['1InsT'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['2DelR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['3DelR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['4DelR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['5DelR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['2InsR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['3InsR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['3InsR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['4InsR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['5InsR'] = [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['2DelM'] = [[0,0]]
				mutations_ID96[sample]['3DelM'] = [[0,0],[0,0]]
				mutations_ID96[sample]['4DelM'] = [[0,0],[0,0],[0,0]]
				mutations_ID96[sample]['5DelM'] = [[0,0],[0,0],[0,0],[0,0],[0,0]]

			for lines in f:
				if pcawg:
					line = lines.strip().split(",")
					line = [x.replace('"','') for x in line]
					if line[1] == 'repeats':
						mut_type = line[2][0] + line[0][0] + line[0][1].lower() + line[0][2].lower() + "R"
					else:
						mut_type = line[2][0] + line[0][0] + line[0][1].lower() + line[0][2].lower() + line[1][0]
					try:
						repeat_size = int(line[3])
					except:
						repeat_size = int(line[3][0])
					if line[1] == 'MH':
						repeat_size -= 1
					sample_index = 4
				else:
					line = lines.strip().split()
					if line[0] not in indel_types_tsb:
						continue
					categories = line[0].split(":")
					bias = categories[0]
					if bias == 'B' or bias == 'N' or bias == 'Q':
						continue
					mut_type = categories[1] + categories[2] + categories[3]

					repeat_size = int(categories[4])
					if categories[3] == 'M':
						repeat_size -= 1
					sample_index = 1

				for sample in samples:
					if mut_type in mutations_ID96[sample].keys():
						if percentage:
							mutCount = float(line[sample_index])
							if mutCount < 1 and mutCount > 0:
								sig_probs = True
						else:
							mutCount = int(line[sample_index])
						if bias == 'T':
							mutations_ID96[sample][mut_type][repeat_size][0] = mutCount
						else:
							mutations_ID96[sample][mut_type][repeat_size][1] = mutCount
					else:
						continue
					sample_index += 1
	except:
		ID96 = False
		print("No ID96 provided")		


########################### ID-simple ################################################
	try:
		indel_types = ['1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5', '1:Del:C:6'
					   '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5', '1:Del:T:6'
					   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
					   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5', 
					   'long_Del', 'long_Ins', 'MH', 'complex']

		with open (sample_matrices_path + "ID/" + project + ".ID28." + file_extension) as f:
			first_line = f.readline()
			for sample in samples:
				mutations_simple[sample] = OrderedDict()
				mutations_simple[sample]['1DelC'] = [0,0,0,0,0,0]
				mutations_simple[sample]['1DelT'] = [0,0,0,0,0,0]
				mutations_simple[sample]['1InsC'] = [0,0,0,0,0,0]
				mutations_simple[sample]['1InsT'] = [0,0,0,0,0,0]
				mutations_simple[sample]['long_Del'] = [0]
				mutations_simple[sample]['long_Ins'] = [0]
				mutations_simple[sample]['MH'] = [0]
				mutations_simple[sample]['complex'] = [0]

			for lines in f:
				line = lines.strip().split()
				categories = line[0].split(":")
				if len(categories) < 2:
					mut_type = categories[0]
					repeat_size = 0
				else:
					mut_type = categories[0] + categories[1] + categories[2]
					repeat_size = int(categories[3])
				sample_index = 1

				for sample in samples:
					#if mut_type in mutations[sample].keys():
						if percentage:
							mutCount = float(line[sample_index])
							if mutCount < 1 and mutCount > 0:
								sig_probs = True
						else:
							mutCount = int(line[sample_index])
						mutations_simple[sample][mut_type][repeat_size] = mutCount

						sample_index += 1
	except:
		ID28 = False
		print("No ID28 provided")



########################### DBS-78 ################################################
	try:
		dinucs = ['TT>GG','TT>CG','TT>AG','TT>GC','TT>CC','TT>AC','TT>GA','TT>CA','TT>AA','AC>CA','AC>CG','AC>CT','AC>GA',
				  'AC>GG','AC>GT','AC>TA','AC>TG','AC>TT','CT>AA','CT>AC','CT>AG','CT>GA','CT>GC','CT>GG','CT>TG','CT>TC',
				  'CT>TA','AT>CA','AT>CC','AT>CG','AT>GA','AT>GC','AT>TA','TG>GT','TG>CT','TG>AT','TG>GC','TG>CC','TG>AC',
				  'TG>GA','TG>CA','TG>AA','CC>AA','CC>AG','CC>AT','CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT','CG>AT',
				  'CG>GC','CG>GT','CG>TC','CG>TA','CG>TT','TC>GT','TC>CT','TC>AT','TC>GG','TC>CG','TC>AG','TC>GA','TC>CA',
				  'TC>AA','GC>AA','GC>AG','GC>AT','GC>CA','GC>CG','GC>TA','TA>GT','TA>CT','TA>AT','TA>GG','TA>CG','TA>GC']

		revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
		mutations = OrderedDict()
		with open (sample_matrices_path + "DBS/" + project + ".DBS78." + file_extension) as f:
			first_line = f.readline()
			for sample in samples:
				mutations_78[sample] = OrderedDict()
				mutations_78[sample]['AC'] = OrderedDict()
				mutations_78[sample]['AT'] = OrderedDict()
				mutations_78[sample]['CC'] = OrderedDict()
				mutations_78[sample]['CG'] = OrderedDict()
				mutations_78[sample]['CT'] = OrderedDict()
				mutations_78[sample]['GC'] = OrderedDict()
				mutations_78[sample]['TA'] = OrderedDict()
				mutations_78[sample]['TC'] = OrderedDict()
				mutations_78[sample]['TG'] = OrderedDict()
				mutations_78[sample]['TT'] = OrderedDict()


			for lines in f:
				if pcawg:
					line = lines.strip().split(",")
					line = [x.replace('"','') for x in line]
					mut = line[0] + ">" + line[1]
					nuc = line[1]
					mut_type = line[0]
					if mut not in dinucs:
						nuc = revcompl(nuc)
						mut_type = revcompl(mut_type)
					sample_index = 2
				else:
					line = lines.strip().split()
					mut = line[0]
					nuc = line[0][3:]
					mut_type = line[0][0:2]
					if mut not in dinucs:
						nuc = revcompl(nuc)
						mut_type = revcompl(mut_type)
					sample_index = 1

				for sample in samples:
					if percentage:
						mutCount = float(line[sample_index])
						if mutCount < 1 and mutCount > 0:
							sig_probs = True
					else:
						mutCount = int(line[sample_index])
					mutations_78[sample][mut_type][nuc] = mutCount
					sample_index += 1
	except:
		DBS78 = False
		print("No DBS78 provided")



########################### DBS-312 ################################################
	try:
		revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
		with open (sample_matrices_path + "DBS/" + project + ".DBS186." + file_extension) as f:
			first_line = f.readline()
			#samples = first_line.strip().split("\t")
			#samples = samples[1:]
			for sample in samples:
				mutations_312[sample] = OrderedDict()
				mutations_312[sample]['CC'] = OrderedDict()
				mutations_312[sample]['CT'] = OrderedDict()
				mutations_312[sample]['TC'] = OrderedDict()
				mutations_312[sample]['TT'] = OrderedDict()

			for lines in f:
				line = lines.strip().split()
				mut = line[0][2:]
				nuc = line[0][5:]
				mut_type = line[0][2:4]
				bias = line[0][0]
				if bias == 'N' or bias == 'B' or bias == 'Q':
					continue
				else:
					if mut not in dinucs:
						if revcompl(mut) not in dinucs:
							continue
						nuc = revcompl(nuc)
						mut_type = revcompl(mut_type)
					sample_index = 1

					for sample in samples:
						if percentage:
							mutCount = float(line[sample_index])
							if mutCount < 1 and mutCount > 0:
								sig_probs = True
						else:
							mutCount = int(line[sample_index])
						if nuc not in mutations_312[sample][mut_type]:
							mutations_312[sample][mut_type][nuc] = [0,0]
						if bias == 'T':
							mutations_312[sample][mut_type][nuc][0] = mutCount
						else:
							mutations_312[sample][mut_type][nuc][1] = mutCount
						sample_index += 1
	except:
		DBS312 = False
		print("No DBS312 provided")	

	########### plot through each sample ####################################################
	for sample in samples:
		plt.rcParams['axes.linewidth'] = 2
		plot1 = plt.figure(figsize=(13.2,8.35))
		plt.rc('axes', edgecolor='lightgray', linewidth=0.5)
		panel1 = plt.axes([0.044, 0.707, 0.172, 0.174]) # SBS-6
		panel2 = plt.axes([0.262, 0.707, 0.144, 0.174]) # SBS-24
		panel3 = plt.axes([0.035, 0.4595, 0.394, 0.097]) # SBS-96
		panel5 = plt.axes([0.035, 0.25, 0.394, 0.097]) # SBS-1536

		panel4 = plt.axes([0.035, 0.05, 0.394, 0.114]) # SBS-384

		panel7 = plt.axes([0.528, 0.776, 0.455, 0.125]) # DBS-78
		#panel7 = plt.axes([0.528, 0.756, 0.455, 0.125]) # DBS-78
		#panel8 = plt.axes([0.783, 0.53, 0.2, 0.15]) # DBS-312
		panel8 = plt.axes([0.783, 0.55, 0.2, 0.15]) # DBS-312

		panel9 = plt.axes([0.528, 0.3, 0.455, 0.155]) # ID-78
		panel10 = plt.axes([0.528, 0.05, 0.455, 0.155]) # ID-96
		# panel11 = plt.axes([0.528, 0.53, 0.15, 0.15]) # ID-simple
		panel11 = plt.axes([0.528, 0.55, 0.15, 0.15]) # ID-simple

		panel12 = plt.axes([0.035, 0.352, 0.394, 0.0485]) # 1536-middle
		panel13 = plt.axes([0.035, 0.406, 0.394, 0.0485]) # 1536-top

		plt.text(.02, .96, sample, fontsize=20, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		# plt.text(.03, .93, "SBS", fontsize=15, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		# plt.text(.5, .94, "DBS", fontsize=15, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		# plt.text(.5, .69, "Indels", fontsize=15, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

		############### plot SBS-96 ##############################################
		if SBS96:
			total_count = sum(sum(nuc.values()) for nuc in mutations_96[sample].values())
			plt.rcParams['axes.linewidth'] = 2
		
			xlabels = []
			
			x = 0.4
			ymax = 0
			colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
			i = 0
			for key in mutations_96[sample]:
				for seq in mutations_96[sample][key]:
					xlabels.append(seq[0]+seq[2]+seq[6])
					if percentage:
						if total_count > 0:	
							panel3.bar(x, mutations_96[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
							if mutations_96[sample][key][seq]/total_count*100 > ymax:
								ymax = mutations_96[sample][key][seq]/total_count*100
					else:
						panel3.bar(x, mutations_96[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
						if mutations_96[sample][key][seq] > ymax:
								ymax = mutations_96[sample][key][seq]
					x += 1
				i += 1

			x = .036
			y3 = .56
			y = int(ymax*1.25)
			y2 = y+2
			for i in range(0, 6, 1):
				panel3.add_patch(plt.Rectangle((x,y3), .064, .005, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
				x += .0657

			yText = y3 + 0.01
			plt.text(.055, yText, 'C>A', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.125, yText, 'C>G', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.19, yText, 'C>T', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.255, yText, 'T>A', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.32, yText, 'T>C', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
			plt.text(.385, yText, 'T>G', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)


			while y%4 != 0:
				y += 1
			ytick_offest = int(y/4)

			if percentage:
				ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
				ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
						  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
			else:
				ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
				ylabels= [0, ytick_offest, ytick_offest*2, 
						  ytick_offest*3, ytick_offest*4]		

			labs = np.arange(0.375,96.375,1)
			if not percentage:
				ylabels = ['{:,}'.format(int(x)) for x in ylabels]
				if len(ylabels[-1]) > 3:
					ylabels_temp = []
					if len(ylabels[-1]) > 7:
						for label in ylabels:
							if len(label) > 7:
								ylabels_temp.append(label[0:-8] + "m")
							elif len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)

					else:
						for label in ylabels:
							if len(label) > 3:
								ylabels_temp.append(label[0:-4] + "k")
							else:
								ylabels_temp.append(label)
					ylabels = ylabels_temp


			# if not percentage:
			# 	ylabels = ['{:,}'.format(int(x)) for x in ylabels]

			panel3.set_xlim([0, 96])
			panel3.set_ylim([0, y])
			panel3.set_xticks(labs)
			panel3.set_yticks(ylabs)
			count = 0
			m = 0

			panel3.set_yticklabels(ylabels, fontsize=8, color='black')
			panel3.yaxis.grid(True)
			panel3.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
			panel3.set_xlabel('')

			if percentage:
				panel3.set_ylabel("Percentage of Single Base Substitutions", fontsize=6, fontname="Times New Roman", weight = 'bold')
			else:
				panel3.set_ylabel("Number of Single Base Substitutions", fontsize=6, fontname="Times New Roman", weight = 'bold')


			panel3.tick_params(axis='both',which='both',\
							   bottom=False, labelbottom=False,\
							   left=True, labelleft=True,\
							   right=True, labelright=False,\
							   top=False, labeltop=False,\
							   direction='in', length=2.5, colors='lightgray', width=0.25)


			[i.set_color("black") for i in panel3.get_yticklabels()]
	
		########################## SBS-6 ###################################################
		total_count = sum(mutations_6[sample].values())
		xlabels = []
		y = -0.5
		xmax = 0
		colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
		i = 0
		for key in mutations_6[sample]:
			xlabels.append(key)
			if percentage:	
				if total_count > 0:
					panel1.barh(y, mutations_6[sample][key]/total_count*100,height=0.7,color=colors[i],align='center', zorder=1000)
					if mutations_6[sample][key]/total_count*100 > xmax:
						xmax = mutations_6[sample][key]/total_count*100
			else:
				panel1.barh(y, mutations_6[sample][key],height=0.7,color=colors[i],align='center', zorder=1000)
				if mutations_6[sample][key] > xmax:
						xmax = mutations_6[sample][key]
			y -= 1
			i += 1

		y = .043
		x3 = .87
		x = int(xmax*1.1)

		while x%4 != 0:
			x += 1
		xtick_offest = int(x/4)

		if percentage:
			xlabs = [0, round(xtick_offest, 1), round(xtick_offest*2, 1), round(xtick_offest*3, 1), round(xtick_offest*4, 1)]
			xlabels= [str(0), str(round(xtick_offest, 1)) + "%", str(round(xtick_offest*2, 1)) + "%", 
					  str(round(xtick_offest*3, 1)) + "%", str(round(xtick_offest*4, 1)) + "%"]
		else:
			xlabs = [0, xtick_offest, xtick_offest*2, xtick_offest*3, xtick_offest*4]
			xlabels= [0, xtick_offest, xtick_offest*2, 
					  xtick_offest*3, xtick_offest*4]	


		if not percentage:
			xlabels = ['{:,}'.format(int(x)) for x in xlabels]
			if len(xlabels[-1]) > 3:
				xlabels_temp = []
				if len(xlabels[-1]) > 7:
					for label in xlabels:
						if len(label) > 7:
							xlabels_temp.append(label[0:-8] + "m")
						elif len(label) > 3:
							xlabels_temp.append(label[0:-4] + "k")
						else:
							xlabels_temp.append(label)

				else:
					for label in xlabels:
						if len(label) > 3:
							xlabels_temp.append(label[0:-4] + "k")
						else:
							xlabels_temp.append(label)
				xlabels = xlabels_temp


		# if not percentage:
		# 	xlabels = ['{:,}'.format(int(x)) for x in xlabels]

		ylabs = np.arange(-5.5, 0.5, 1)
		ylabels = (['T>G','T>C','T>A','C>T','C>G','C>A'])
		panel1.set_xlim([0, x])
		panel1.set_ylim([-6, 0])
		panel1.set_xticks(xlabs)
		panel1.set_yticks(ylabs)
		panel1.set_xticklabels(xlabels, fontsize=8)
		panel1.set_yticklabels(ylabels, fontsize=8)
		panel1.spines['right'].set_visible(False)
		panel1.spines['top'].set_visible(False)

		# if sig_probs:
		panel1.text(.05, .9, "SBS-6", fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		# else:
		# 	panel1.text(.045, .9, sample + ": " + "{:,}".format(int(total_count)) + " subs", fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)


		if percentage:
			panel1.set_xlabel("Percentage of Single Base Substitutions", fontsize=10, fontname="Times New Roman", weight = 'bold')
		else:
			panel1.set_xlabel("Number of Single Base Substitutions", fontsize=10, fontname="Times New Roman", weight = 'bold')

		panel1.set_ylabel('')

		panel1.tick_params(axis='both',which='both',\
						   bottom=True, labelbottom=True,\
						   left=False, labelleft=True,\
						   right=False, labelright=False,\
						   top=False, labeltop=False,\
						   width=0.2, color='lightgrey')


		########################## SBS-24 ###################################################
		total_count = sum(sum(tsb) for tsb in mutations_24[sample].values())
		y = 12.485
		xmax = 0
		for key in mutations_24[sample]:
			if percentage:
				if total_count > 0:
					trans = panel2.barh(y, mutations_24[sample][key][0]/total_count*100,height=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
					y -= 0.75
					untrans = panel2.barh(y, mutations_24[sample][key][1]/total_count*100,height=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
					y -= .2475
					if mutations_24[sample][key][0]/total_count*100 > xmax:
							xmax = mutations_24[sample][key][0]/total_count*100
					if mutations_24[sample][key][1]/total_count*100 > xmax:
							xmax = mutations_24[sample][key][1]/total_count*100

			else:
				trans = panel2.barh(y, mutations_24[sample][key][0],height=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
				y -= 0.75
				untrans = panel2.barh(y, mutations_24[sample][key][1],height=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
				y -= .2475
				if mutations_24[sample][key][0] > xmax:
						xmax = mutations_24[sample][key][0]
				if mutations_24[sample][key][1] > xmax:
						xmax = mutations_24[sample][key][1]
			y -= 1

		x = int(xmax*1.1)

		while x%4 != 0:
			x += 1

		xtick_offest = int(x/4)

			
		if percentage:
			xlabs = [0, round(xtick_offest, 1), round(xtick_offest*2, 1), round(xtick_offest*3, 1), round(xtick_offest*4, 1)]
			xlabels= [str(0), str(round(xtick_offest, 1)) + "%", str(round(xtick_offest*2, 1)) + "%", 
					  str(round(xtick_offest*3, 1)) + "%", str(round(xtick_offest*4, 1)) + "%"]
		else:
			xlabs = [0, xtick_offest, xtick_offest*2, xtick_offest*3, xtick_offest*4]
			xlabels= [0, xtick_offest, xtick_offest*2, 
					  xtick_offest*3, xtick_offest*4]
		if not percentage:
			xlabels = ['{:,}'.format(int(x)) for x in xlabels]
			if len(xlabels[-1]) > 3:
				xlabels_temp = []
				if len(xlabels[-1]) > 7:
					for label in xlabels:
						if len(label) > 7:
							xlabels_temp.append(label[0:-8] + "m")
						elif len(label) > 3:
							xlabels_temp.append(label[0:-4] + "k")
						else:
							xlabels_temp.append(label)

				else:
					for label in xlabels:
						if len(label) > 3:
							xlabels_temp.append(label[0:-4] + "k")
						else:
							xlabels_temp.append(label)
				xlabels = xlabels_temp


		# if not percentage:
		# 	xlabels = ['{:,}'.format(int(x)) for x in xlabels]
		ylabs = np.arange(2.15, 13, 2)
		ylabels = (['T>G','T>C','T>A','C>T','C>G','C>A'])
		panel2.set_xlim([0, x])
		panel2.set_ylim([1.2524, 13.235])
		panel2.set_yticks(ylabs)
		panel2.set_xticks(xlabs)
		panel2.set_xticklabels(xlabels, fontsize=10)
		panel2.set_yticklabels(ylabels, fontsize=10)
		panel2.set_xlabel('')
		panel2.set_ylabel('')
		panel2.spines['right'].set_visible(False)
		panel2.spines['top'].set_visible(False)

		# if sig_probs:
		panel2.text(.27, .9, "SBS-24", fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		# else:
		# 	panel2.text(.265, .9, sample + ": " + "{:,}".format(int(total_count)) + " transcribed subs", fontsize=8, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

		if percentage:
			panel2.set_xlabel("Percentage of Single Base Substitutions", fontsize=10, fontname="Times New Roman", weight = 'bold')
		else:
			panel2.set_xlabel("Number of Single Base Substitutions", fontsize=10, fontname="Times New Roman", weight = 'bold')

		panel2.tick_params(axis='both',which='both',\
						   bottom=True, labelbottom=True,\
						   left=False, labelleft=True,\
						   right=False, labelright=False,\
						   top=False, labeltop=False,\
						   width=0.2, color='lightgrey')
		# panel2.legend(handles=[trans, untrans], prop={'size':25})

		########################## SBS-384 ###################################################
		total_count = sum(sum(sum(tsb) for tsb in nuc.values()) for nuc in mutations_384[sample].values())
		xlabels = []
		
		x = 0.7
		ymax = 0
		colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
		i = 0
		for key in mutations_384[sample]:
			for seq in mutations_384[sample][key]:
				xlabels.append(seq[0]+seq[2]+seq[6])
				if percentage:
					if total_count > 0:
						trans = panel4.bar(x, mutations_384[sample][key][seq][0]/total_count*100,width=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
						x += 0.75
						untrans = panel4.bar(x, mutations_384[sample][key][seq][1]/total_count*100,width=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
						x += .2475
						if mutations_384[sample][key][seq][0]/total_count*100 > ymax:
								ymax = mutations_384[sample][key][seq][0]/total_count*100
						if mutations_384[sample][key][seq][1]/total_count*100 > ymax:
								ymax = mutations_384[sample][key][seq][1]/total_count*100

				else:
					trans = panel4.bar(x, mutations_384[sample][key][seq][0],width=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
					x += 0.75
					untrans = panel4.bar(x, mutations_384[sample][key][seq][1],width=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
					x += .2475
					if mutations_384[sample][key][seq][0] > ymax:
							ymax = mutations_384[sample][key][seq][0]
					if mutations_384[sample][key][seq][1] > ymax:
							ymax = mutations_384[sample][key][seq][1]
				x += 1
			i += 1

		x = .036
		# y3 = .2775

		y3 = .1675
		y = int(ymax*1.25)
		x_plot = 0
		while y%4 != 0:
			y += 1


		yText = y3 + 0.01
		plt.text(.055, yText, 'C>A', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.125, yText, 'C>G', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.19, yText, 'C>T', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.255, yText, 'T>A', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.32, yText, 'T>C', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.385, yText, 'T>G', fontsize=10, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

		ytick_offest = int(y/4)

		for i in range(0, 6, 1):
			panel4.add_patch(plt.Rectangle((x,y3), .064, .005, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			panel4.add_patch(plt.Rectangle((x_plot,0), 32, y, facecolor=colors[i], zorder=0, alpha = 0.25, edgecolor='grey'))
			x += .0657
			x_plot += 32
			
		if percentage:
			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
		else:
			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
			ylabels= [0, ytick_offest, ytick_offest*2, 
					  ytick_offest*3, ytick_offest*4]

		labs = np.arange(0.750,192.750,1)
		if not percentage:
			ylabels = ['{:,}'.format(int(x)) for x in ylabels]
			if len(ylabels[-1]) > 3:
				ylabels_temp = []
				if len(ylabels[-1]) > 7:
					for label in ylabels:
						if len(label) > 7:
							ylabels_temp.append(label[0:-8] + "m")
						elif len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)

				else:
					for label in ylabels:
						if len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)
				ylabels = ylabels_temp
		# if not percentage:
		# 	ylabels = ['{:,}'.format(int(x)) for x in ylabels]
		panel4.set_xlim([0, 96])
		panel4.set_ylim([0, y])
		panel4.set_xticks(labs)
		panel4.set_yticks(ylabs)
		count = 0
		m = 0
		for i in range (0, 96, 1):
			panel4.text(i/243.5 + .0355, .035, xlabels[i][0], fontsize=4, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			panel4.text(i/243.5 + .0355, .04, xlabels[i][1], fontsize=4, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
			panel4.text(i/243.5 + .0355, .045, xlabels[i][2], fontsize=4, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			count += 1
			if count == 16:
				count = 0
				m += 1
		
		plt.text(0.04, 0.2, "SBS-384", fontsize=10, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

		panel4.set_yticklabels(ylabels, fontsize=8)
		panel4.yaxis.grid(True)
		panel4.grid(which='major', axis='y', color=[0.7,0.7,0.7], zorder=1)
		panel4.legend(handles=[trans, untrans], prop={'size':4})
		if percentage:
			panel4.set_ylabel("Percentage of Single Base Substitutions", fontsize=6, fontname="Times New Roman", weight = 'bold')
		else:
			panel4.set_ylabel("Number of Single Base Substitutions", fontsize=6, fontname="Times New Roman", weight = 'bold')

		panel4.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=False,\
						   left=False, labelleft=True,\
						   right=False, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors=[0.6, 0.6, 0.6])

		[i.set_color("black") for i in panel4.get_yticklabels()]




	############### plot SBS-1536 ##############################################
		xlabels = []
		ylabels = []
		ylabels_5 = []
		ylabels_3 = []
		

		colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
		colors_heat = [np.linspace(56/255,255/255,5), np.linspace(66/255,225/255,5), np.linspace(157/255,40/255,5)]
		colors_heat_compact = [np.linspace(56/255,255/255,5), np.linspace(66/255,225/255,5), np.linspace(157/255,40/255,5)]

		total_count_sample = sum(sum(nuc.values()) for nuc in mutations_96[sample].values())
		total_count = max_all[sample]*1.1
		ratio = total_count/total_count_sample

		i = 0
		x_pos = 0
		x_inter = 0
		for key in mutations_1536[sample]:
			y_pos = 15
			for penta in mutations_1536[sample][key]:
				key_5 = penta[0]
				key_3 = penta[1]
				ylabels.append(penta[0] + "---" + penta[1])
				for tri in mutations_1536[sample][key][penta]:
					tri_nuc = tri[0] + "[" + key + "]" + tri[1]
					normalized = mutations_96[sample][key][tri_nuc]
					try:
						mut_count = int(int(20 * round(float(mutations_1536[sample][key][penta][tri]/total_count_sample/ratio* 100))/20)/20)
						mutations_5[sample][key][key_5][tri] += float(mutations_1536[sample][key][penta][tri])
						mutations_3[sample][key][key_3][tri] += float(mutations_1536[sample][key][penta][tri])
						if mutations_5[sample][key][key_5][tri] > max_5[sample]:
							max_5[sample] = mutations_5[sample][key][key_5][tri]
						if mutations_3[sample][key][key_3][tri] > max_3[sample]:
							max_3[sample] = mutations_3[sample][key][key_3][tri]

					except:
						mut_count = 0
					xlabels.append(tri[0]+"-"+tri[1])
					rectangle=mplpatches.Rectangle((x_pos, y_pos), 1, 1,\
													linewidth=1,\
													facecolor=(colors_heat[0][mut_count], colors_heat[1][mut_count], colors_heat[2][mut_count]))

					panel5.add_patch(rectangle)
					x_pos += 1
				y_pos -= 1
				x_pos = x_inter

			x_inter += 17
			x_pos = x_inter
			i += 1
	


		x_pos = 0 
		x_inter = 0
		total_count_5 = max_5[sample]*1.1
		total_count_3 = max_3[sample]*1.1
		ratio_5 = total_count_5/total_count_sample
		ratio_3 = total_count_3/total_count_sample
		ratio_total = max(ratio_5, ratio_3)

		for key in mutations_5[sample]:
			y_pos = 3
			for penta in mutations_5[sample][key]:
				ylabels_5.append(penta + "---N" )
				ylabels_3.append("N---" + penta)
				for tri in mutations_5[sample][key][penta]:
					mut_count = int(int(20 * round(float(mutations_5[sample][key][penta][tri])/total_count_sample/ratio_total*100)/20)/20)
					mut_count_3 = int(int(20 * round(float(mutations_3[sample][key][penta][tri])/total_count_sample/ratio_total*100)/20)/20)
					rectangle=mplpatches.Rectangle((x_pos, y_pos), 1, 1,\
													linewidth=1,\
													facecolor=(colors_heat_compact[0][mut_count], colors_heat_compact[1][mut_count], colors_heat_compact[2][mut_count]))

					panel13.add_patch(rectangle)
					rectangle=mplpatches.Rectangle((x_pos, y_pos), 1, 1,\
													linewidth=1,\
													facecolor=(colors_heat[0][mut_count_3], colors_heat[1][mut_count_3], colors_heat[2][mut_count_3]))

					panel12.add_patch(rectangle)
					x_pos += 1
				y_pos -= 1
				x_pos = x_inter
			x_inter += 17
			x_pos = x_inter
			i += 1

		x = 0.5
		ymax = 0
		colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
		i = 0


		# scale bar for bottom 1536 plot
		y_grad = .0955/len(colors_heat[0])
		y_start = 0.251
		for l in range (0, len(colors_heat[0]), 1):
			rectangle=mplpatches.Rectangle((.4325, y_start), .005, y_grad,\
											linewidth=1,\
											facecolor=(colors_heat[0][l], colors_heat[1][l], colors_heat[2][l]),\
											transform=plt.gcf().transFigure, clip_on=False,\
											edgecolor=(colors_heat[0][l], colors_heat[1][l], colors_heat[2][l]))
			panel5.add_patch(rectangle)
			y_start += y_grad

		# scale bar for top 1536 plot
		y_grad = .0474/len(colors_heat_compact[0])
		y_start = 0.406
		for l in range (0, len(colors_heat_compact[0]), 1):
			rectangle=mplpatches.Rectangle((.4325, y_start), .005, y_grad,\
											linewidth=1,\
											facecolor=(colors_heat_compact[0][l], colors_heat_compact[1][l], colors_heat_compact[2][l]),\
											transform=plt.gcf().transFigure, clip_on=False,\
											edgecolor=(colors_heat_compact[0][l], colors_heat_compact[1][l], colors_heat_compact[2][l]))
			panel1.add_patch(rectangle)
			y_start += y_grad

		# scale bar for middle 1536 plot
		y_grad = .0474/len(colors_heat_compact[0])
		y_start = 0.3525
		for l in range (0, len(colors_heat_compact[0]), 1):
			rectangle=mplpatches.Rectangle((.4325, y_start), .005, y_grad,\
											linewidth=1,\
											facecolor=(colors_heat_compact[0][l], colors_heat_compact[1][l], colors_heat_compact[2][l]),\
											transform=plt.gcf().transFigure, clip_on=False,\
											edgecolor=(colors_heat_compact[0][l], colors_heat_compact[1][l], colors_heat_compact[2][l]))
			panel1.add_patch(rectangle)
			y_start += y_grad
		# y_tick_grad = max_count[sample]/2
	
		# scale numbers for bottom 1536 plot
		plt.text(.439, .251, '0', fontsize=5, fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.439, .3, str(ratio/2)[:5], fontsize=5, fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.439, .34, str(ratio)[:5], fontsize=5, fontweight='bold', transform=plt.gcf().transFigure)
		y = int(ymax*1.25)

		# scale numbers for top 1536 plot
		plt.text(.439, .406, '0', fontsize=5, fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.439, .426, str(ratio_total/2)[:5], fontsize=5, fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.439, .448, str(ratio_total)[:5], fontsize=5, fontweight='bold', transform=plt.gcf().transFigure)
		y = int(ymax*1.25)

		# scale numbers for middle 1536 plot
		plt.text(.439, .3525, '0', fontsize=5, fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.439, .375, str(ratio_total/2)[:5], fontsize=5, fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.439, .395, str(ratio_total)[:5], fontsize=5, fontweight='bold', transform=plt.gcf().transFigure)
		y = int(ymax*1.25)

		if y <= 4:
			y += 4

		while y%4 != 0:
			y += 1
		ytick_offest = int(y/4)


		panel5.set_ylim([0, 16])
		panel5.set_xlim([0, 101])
		panel13.set_xlim([0, 101])
		panel13.set_ylim([0, 4])
		panel12.set_xlim([0, 101])
		panel12.set_ylim([0, 4])
		panel12.set_yticks([])
		panel13.set_yticks([])
		panel13.set_xticks([])
		panel13.set_xticks([])
		panel5.set_xticks([])
		panel5.set_yticks([])


		# x-axis 1536 bottom plot
		m = 0
		count = 0
		x_letter = 0
		for i in range (0, 96, 1):
			plt.text(x_letter/236 + .0355, .237, xlabels[i][0], fontsize=4, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			plt.text(x_letter/236 + .0355, .241, xlabels[i][1], fontsize=4, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
			plt.text(x_letter/236+ .0355, .245, xlabels[i][2], fontsize=4, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			count += 1
			x_letter += .92
			if (i+1)%16 == 0 and i != 0:
				x_letter += .92
			if count == 16:
				count = 0
				m += 1	

		# y-axis 1536 bottom plot
		m = 0
		count = 0
		y_letter = 5.5
		for i in range (0, 16, 1):
			plt.text(.019, y_letter/16 , ylabels[i][0], fontsize=4, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			plt.text(.022, y_letter/16, ylabels[i][1:4], fontsize=4, color='black', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
			plt.text(.03, y_letter/16, ylabels[i][4], fontsize=4, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			count += 1
			y_letter -= .0975
			
			if count == 16:
				count = 0
				m += 1	

		# y-axis 1536 top matrix plot
		m = 0
		count = 0
		y_letter = 7.17
		for i in range (0, 4, 1):
			plt.text(.019, y_letter/16 , ylabels_5[i][0], fontsize=4, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			plt.text(.022, y_letter/16 + 0, ylabels_5[i][1:4], fontsize=4, color='black', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
			plt.text(.03, y_letter/16 + 0, ylabels_5[i][4], fontsize=4, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			count += 1
			y_letter -= .195
			
			if count == 16:
				count = 0
				m += 1	

		# y-axis 1536 middle matrix plot
		m = 0
		count = 0
		y_letter = 6.3
		for i in range (0, 4, 1):
			plt.text(.019, y_letter/16 , ylabels_3[i][0], fontsize=4, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			plt.text(.022, y_letter/16 + 0, ylabels_3[i][1:4], fontsize=4, color='black', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
			plt.text(.03, y_letter/16 + 0, ylabels_3[i][4], fontsize=4, color='black', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			count += 1
			y_letter -= .195
			
			if count == 16:
				count = 0
				m += 1	




		panel2.text(0.04, 0.6, "SBS-96 and SBS-1536", fontsize=10, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)



		panel5.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=False,\
						   left=False, labelleft=False,\
						   right=True, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors='white', width=2)
		panel13.axis('off')
		panel13.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=False,\
						   left=False, labelleft=False,\
						   right=True, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors='white', width=2)
		panel12.axis('off')
		panel12.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=False,\
						   left=False, labelleft=False,\
						   right=True, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors='white', width=2)

		[i.set_color("black") for i in plt.gca().get_yticklabels()]




	############### plot ID-83 ##############################################

		xlabels = []
		
		x = 0.4
		ymax = 0
		colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256], 
				  [253/256,202/256,181/256], [252/256,138/256,106/256], [241/256,68/256,50/256], [188/256,25/256,26/256],
				  [208/256,225/256,242/256], [148/256,196/256,223/256], [74/256,152/256,201/256], [23/256,100/256,171/256],
				  [226/256,226/256,239/256], [182/256,182/256,216/256], [134/256,131/256,189/256], [98/256,64/256,155/256]]

		i = 0
		for key in mutations_83[sample]:
			l = 1
			for seq in mutations_83[sample][key]:
				xlabels.append(l)
				if percentage:
					if total_count > 0:
						panel9.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
						if seq/total_count*100 > ymax:
							ymax = seq/total_count*100
				else:
					panel9.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq > ymax:
							ymax = seq
				x += 1
				l += 1
			i += 1

		x = .5292
		y_top = .4575
		y_bottom = .2875
		y = int(ymax*1.25)
		y2 = y+2
		for i in range(0, 12, 1):
			panel9.add_patch(plt.Rectangle((x,y_top), .02935, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			panel9.add_patch(plt.Rectangle((x,y_bottom), .02935, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			x += .0329

		panel9.add_patch(plt.Rectangle((x-.0005,y_top), .003, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
		panel9.add_patch(plt.Rectangle((x-.0005,y_bottom), .003, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
		x +=.00525
		panel9.add_patch(plt.Rectangle((x,y_top), .00775, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
		panel9.add_patch(plt.Rectangle((x,y_bottom), .00775, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
		x += .011
		panel9.add_patch(plt.Rectangle((x,y_top), .01315, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
		panel9.add_patch(plt.Rectangle((x,y_bottom), .01315, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))		
		x += .0165
		panel9.add_patch(plt.Rectangle((x,y_top), .024, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
		panel9.add_patch(plt.Rectangle((x,y_bottom), .024, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))




		yText = y_top + 0.001
		plt.text(.5425, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.5745, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
		plt.text(.6065, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.64, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
		plt.text(.6735, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.7065, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.739, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.771, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		plt.text(.805, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.838, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.871, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.903, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		plt.text(.9235, yText, '2', fontsize=6, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.9313, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.945, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.965, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

		yText_labels_top = yText + 0.02
		yText_labels_bottom = y_bottom - .01
		yText_labels_bottom_sec = yText_labels_bottom - .01

		plt.text(.5425, yText_labels_top, '1bp Deletion', fontsize=7, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
		plt.text(.6065, yText_labels_top, '1bp Insertion', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.689, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.82, yText_labels_top, '>1bp Insertions at Repeats\n       (Insertion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.926, yText_labels_top, ' Mircohomology\n(Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		plt.text(.535, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=5, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
		plt.text(.6, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.698, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.83, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.925, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x = .5292
		for i in range (0, 8, 1):
			if i != 2 and i != 3:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

			x += .03285

		x += 0.0005
		for i in range (0, 4, 1):
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			x += .033

		x -= 0.0009
		plt.text(x, yText_labels_bottom, '1', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .006
		plt.text(x, yText_labels_bottom, '1  2', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .011
		plt.text(x, yText_labels_bottom, '1  2  3', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .0165
		plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)


		while y%4 != 0:
			y += 1
		ytick_offest = int(y/4)

		if percentage:
			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
		else:
			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
			ylabels= [0, ytick_offest, ytick_offest*2, 
					  ytick_offest*3, ytick_offest*4]

		labs = np.arange(0.375,83.375,1)
		# if not percentage:
		# 	ylabels = ['{:,}'.format(int(x)) for x in ylabels]
		if not percentage:
			ylabels = ['{:,}'.format(int(x)) for x in ylabels]
			if len(ylabels[-1]) > 3:
				ylabels_temp = []
				if len(ylabels[-1]) > 7:
					for label in ylabels:
						if len(label) > 7:
							ylabels_temp.append(label[0:-8] + "m")
						elif len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)

				else:
					for label in ylabels:
						if len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)
				ylabels = ylabels_temp

		panel9.set_xlim([0, 83])
		panel9.set_ylim([0, y])
		panel9.set_xticks(labs)
		panel9.set_yticks(ylabs)	

		plt.text(0.528, 0.48 + 0.01, "ID-83", fontsize=10, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

		panel9.set_yticklabels(ylabels, fontsize=8)
		panel9.yaxis.grid(True)
		panel9.grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
		panel9.set_xlabel('')
		panel9.set_ylabel('')

		if percentage:
			panel9.set_ylabel("Percentage of Indels", fontsize=8, fontname="Times New Roman", weight = 'bold')
		else:
			panel9.set_ylabel("Number of Indels", fontsize=8, fontname="Times New Roman", weight = 'bold')

		panel9.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=False,\
						   left=False, labelleft=True,\
						   right=False, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors='gray', width=2)

		[i.set_color("black") for i in panel9.get_yticklabels()]


	############### plot ID-96 ##############################################
		xlabels = []
		total_count = sum(sum(sum(i) for i in nuc) for nuc in mutations_ID96[sample].values())
		
		x = 0.4
		ymax = 0
		colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256], 
				  [253/256,202/256,181/256], [252/256,138/256,106/256], [241/256,68/256,50/256], [188/256,25/256,26/256],
				  [208/256,225/256,242/256], [148/256,196/256,223/256], [74/256,152/256,201/256], [23/256,100/256,171/256],
				  [226/256,226/256,239/256], [182/256,182/256,216/256], [134/256,131/256,189/256], [98/256,64/256,155/256]]

		i = 0
		for key in mutations_ID96[sample]:
			l = 1
			for seq in mutations_ID96[sample][key]:
				xlabels.append(l)
				if percentage:
					if total_count > 0:
						trans = panel10.bar(x, seq[0]/total_count*100,width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
						x += 0.2
						untrans = panel10.bar(x, seq[1]/total_count*100,width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
						if seq[0]/total_count*100 > ymax:
								ymax = seq[0]/total_count*100
						if seq[1]/total_count*100 > ymax:
								ymax = seq[1]/total_count*100

				else:
					trans = panel10.bar(x, seq[0],width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
					x += 0.2
					untrans = panel10.bar(x, seq[1],width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
					if seq[0] > ymax:
							ymax = seq[0]
					if seq[1] > ymax:
							ymax = seq[1]

				x += 0.799
				l += 1
			i += 1

		x = .529
		y_top = .2075
		y_bottom = .0375
		y = int(ymax*1.25)
		y2 = y+2
		for i in range(0, 12, 1):
			panel10.add_patch(plt.Rectangle((x,y_top), .02935, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			panel10.add_patch(plt.Rectangle((x,y_bottom), .02935, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			x += .0329


		panel10.add_patch(plt.Rectangle((x-.0005,y_top), .003, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
		panel10.add_patch(plt.Rectangle((x-.0005,y_bottom), .003, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
		x +=.00525
		panel10.add_patch(plt.Rectangle((x,y_top), .00775, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
		panel10.add_patch(plt.Rectangle((x,y_bottom), .00775, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
		x += .011
		panel10.add_patch(plt.Rectangle((x,y_top), .01315, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
		panel10.add_patch(plt.Rectangle((x,y_bottom), .01315, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))		
		x += .0165
		panel10.add_patch(plt.Rectangle((x,y_top), .024, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
		panel10.add_patch(plt.Rectangle((x,y_bottom), .024, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))



		yText = y_top + 0.001
		plt.text(.5425, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.5745, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
		plt.text(.6065, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.64, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
		plt.text(.6735, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.7065, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.739, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.771, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		plt.text(.805, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.838, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.871, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.903, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		plt.text(.9235, yText, '2', fontsize=6, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.9313, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.945, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.965, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

		yText_labels_top = yText + 0.02
		yText_labels_bottom = y_bottom - .01
		yText_labels_bottom_sec = yText_labels_bottom - .01

		plt.text(.5425, yText_labels_top, '1bp Deletion', fontsize=7, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
		plt.text(.6065, yText_labels_top, '1bp Insertion', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.689, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.82, yText_labels_top, '>1bp Insertions at Repeats\n       (Insertion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.926, yText_labels_top, ' Mircohomology\n(Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		plt.text(.535, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=5, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
		plt.text(.6, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.698, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.83, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.925, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x = .5292
		for i in range (0, 8, 1):
			if i != 2 and i != 3:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

			x += .03285

		x += 0.0005
		for i in range (0, 4, 1):
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			x += .033

		x -= 0.0009
		plt.text(x, yText_labels_bottom, '1', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .006
		plt.text(x, yText_labels_bottom, '1  2', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .011
		plt.text(x, yText_labels_bottom, '1  2  3', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .0165
		plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)



		if y <= 4:
			y += 4

		while y%4 != 0:
			y += 1
		ytick_offest = int(y/4)

		if percentage:
			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
		else:
			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
			ylabels= [0, ytick_offest, ytick_offest*2, 
					  ytick_offest*3, ytick_offest*4]


		if not percentage:
			ylabels = ['{:,}'.format(int(x)) for x in ylabels]
			if len(ylabels[-1]) > 3:
				ylabels_temp = []
				if len(ylabels[-1]) > 7:
					for label in ylabels:
						if len(label) > 7:
							ylabels_temp.append(label[0:-8] + "m")
						elif len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)

				else:
					for label in ylabels:
						if len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)
				ylabels = ylabels_temp

		panel10.set_xlim([0, 83])
		panel10.set_ylim([0, y])
		panel10.set_xticks(labs)
		panel10.set_yticks(ylabs)	

		plt.text(0.528, 0.24, "ID-415", fontsize=10, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)



		panel10.set_yticklabels(ylabels, fontsize=8, color='b')
		panel10.yaxis.grid(True)
		panel10.grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
		panel10.set_xlabel('')
		panel10.set_ylabel('')
		panel10.legend(handles=[trans, untrans], prop={'size':4})

		if percentage:
			panel10.set_ylabel("Percentage of Indels", fontsize=8, fontname="Times New Roman", weight = 'bold')
		else:
			panel10.set_ylabel("Number of Indels", fontsize=8, fontname="Times New Roman", weight = 'bold')

		panel10.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=False,\
						   left=False, labelleft=True,\
						   right=False, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors='gray', width=2)

		[i.set_color("black") for i in panel10.get_yticklabels()]

	# ############### plot ID-simple ##############################################
		total_count = sum(sum(nuc) for nuc in mutations_simple[sample].values())
		xlabels = []
		
		x = 0.4
		ymax = 0
		colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256], 
				  [23/256,100/256,171/256],[98/256,64/256,155/256], [98/256,64/256,155/256]]

		i = 0
		for key in mutations_simple[sample]:
			l = 1
			for seq in mutations_simple[sample][key]:
				xlabels.append(l)
				if percentage:
					if total_count > 0:
						panel11.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
						if seq/total_count*100 > ymax:
							ymax = seq/total_count*100
				else:
					panel11.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq > ymax:
							ymax = seq
				x += 1
				l += 1
			if i < 4:
				i += 1
		x = .529
		y_top = .6825 + 0.02
		y_bottom = .5175 + 0.02
		y = int(ymax*1.25)
		y2 = y+2
	
		for i in range(0, 4, 1):
			panel11.add_patch(plt.Rectangle((x,y_top), .029, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			panel11.add_patch(plt.Rectangle((x,y_bottom), .029, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			x += .032

		panel1.add_patch(plt.Rectangle((x,y_top), .02, .01, facecolor=colors[i+1], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((x,y_bottom), .02, .01, facecolor=colors[i+1], clip_on=False, transform=plt.gcf().transFigure)) 

		yText = y_top + 0.001
		plt.text(.542, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.573, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
		plt.text(.605, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.6375, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)

		yText_labels_top = yText + .012
		yText_labels_bottom = y_bottom 
		yText_labels_bottom_sec = yText_labels_bottom -.015

		plt.text(.545, yText_labels_top, '1bp Deletion', fontsize=5, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
		plt.text(.61, yText_labels_top, '1bp Insertion', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.535, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=5, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
		plt.text(.6, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.66, yText_labels_top, '>1bp', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.661, yText_labels_bottom_sec, 'Type', fontsize=5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x = .529
		yText_labels_bottom = y_bottom -.006#+0.002

		for l in range (0, 4, 1):
			if l < 2:
				for i in range(1, 6, 1):
					plt.text(x, yText_labels_bottom, str(i), fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
					x += 0.0054
				x -= 0.003
				plt.text(x, yText_labels_bottom, '+6', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				x += 0.0082
			else:
				if l == 2:
					x += 0
				for i in range(0, 5, 1):
					plt.text(x, yText_labels_bottom, str(i), fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
					x += 0.0054
				x -= 0.003
				plt.text(x, yText_labels_bottom, '+5', fontsize=4.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
				x += 0.0082

		yText_labels_bottom += 0.002
		plt.text(x, yText_labels_bottom, 'Del', fontsize=3.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure, rotation='vertical')
		x += 0.005
		plt.text(x, yText_labels_bottom, 'Ins', fontsize=3.8, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure, rotation='vertical')
		x += 0.005
		yText_labels_bottom += 0.0005
		plt.text(x, yText_labels_bottom, 'MH', fontsize=3.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure, rotation='vertical')
		x += 0.006
		yText_labels_bottom += 0.001
		plt.text(x, yText_labels_bottom, 'COMP', fontsize=2, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure, rotation='vertical')
	

		while y%4 != 0:
			y += 1
		ytick_offest = int(y/4)

		if percentage:
			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
		else:
			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
			ylabels= [0, ytick_offest, ytick_offest*2, 
					  ytick_offest*3, ytick_offest*4]
		if not percentage:
			ylabels = ['{:,}'.format(int(x)) for x in ylabels]
			if len(ylabels[-1]) > 3:
				ylabels_temp = []
				if len(ylabels[-1]) > 7:
					for label in ylabels:
						if len(label) > 7:
							ylabels_temp.append(label[0:-8] + "m")
						elif len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)

				else:
					for label in ylabels:
						if len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)
				ylabels = ylabels_temp
		# if not percentage:
		# 	ylabels = ['{:,}'.format(int(x)) for x in ylabels]
		panel11.set_xlim([0, 28])
		panel11.set_ylim([0, y])
		panel11.set_yticks(ylabs)	

		plt.text(0.528, 0.71 + 0.02, "ID-28", fontsize=10, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

		panel11.set_yticklabels(ylabels, fontsize=8)
		panel11.yaxis.grid(True)
		panel11.grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
		panel11.set_xlabel('')
		panel11.set_ylabel('')

		if percentage:
			panel11.set_ylabel("Percentage of Indels", fontsize=8, fontname="Times New Roman", weight = 'bold')
		else:
			panel11.set_ylabel("Number of Indels", fontsize=8, fontname="Times New Roman", weight = 'bold')
		panel11.axis('on')
		panel11.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=False,\
						   left=False, labelleft=True,\
						   right=False, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors='gray', width=2)

		[i.set_color("black") for i in panel11.get_yticklabels()]




	############### plot DBS-78 ##############################################
		total_count = sum(sum(nuc.values()) for nuc in mutations_78[sample].values())
		xlabels = []	
		x = 0.4
		ymax = 0
		colors = [[3/256,189/256,239/256], [3/256,102/256,204/256],[162/256,207/256,99/256], 
				  [1/256,102/256,1/256], [255/256,153/256,153/256], [228/256,41/256,38/256], 
				  [255/256,178/256,102/256], [255/256,128/256,1/256], [204/256,153/256,255/256], 
				  [76/256,1/256,153/256]]
		i = 0
		for key in mutations_78[sample]:
			muts = mutations_78[sample][key].keys()
			muts = sorted(muts)
			for seq in muts:
				xlabels.append(seq)
				if percentage:	
					if total_count > 0:
						panel7.bar(x, mutations_78[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
						if mutations_78[sample][key][seq]/total_count*100 > ymax:
								ymax = mutations_78[sample][key][seq]/total_count*100					
				else:
					panel7.bar(x, mutations_78[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations_78[sample][key][seq] > ymax:
							ymax = mutations_78[sample][key][seq]
				x += 1
			i += 1

		x = .5
		y3 = .8825 + 0.02
		y = ymax*1.25
		y2 = y+2
		i = 0
		panel7.add_patch(plt.Rectangle((.53,y3), .048, .01, facecolor=colors[0], clip_on=False, transform=plt.gcf().transFigure)) 
		panel7.add_patch(plt.Rectangle((.581,y3), .032, .01, facecolor=colors[1], clip_on=False, transform=plt.gcf().transFigure)) 
		panel7.add_patch(plt.Rectangle((.61675,y3), .048, .01, facecolor=colors[2], clip_on=False, transform=plt.gcf().transFigure)) 
		panel7.add_patch(plt.Rectangle((.6685,y3), .032, .01, facecolor=colors[3], clip_on=False, transform=plt.gcf().transFigure)) 
		panel7.add_patch(plt.Rectangle((.704,y3), .049, .01, facecolor=colors[4], clip_on=False, transform=plt.gcf().transFigure)) 
		panel7.add_patch(plt.Rectangle((.757,y3), .031, .01, facecolor=colors[5], clip_on=False, transform=plt.gcf().transFigure)) 
		panel7.add_patch(plt.Rectangle((.792,y3), .031, .01, facecolor=colors[6], clip_on=False, transform=plt.gcf().transFigure)) 
		panel7.add_patch(plt.Rectangle((.827,y3), .049, .01, facecolor=colors[7], clip_on=False, transform=plt.gcf().transFigure)) 
		panel7.add_patch(plt.Rectangle((.88,y3), .048, .01, facecolor=colors[8], clip_on=False, transform=plt.gcf().transFigure)) 
		panel7.add_patch(plt.Rectangle((.932,y3), .048, .01, facecolor=colors[9], clip_on=False, transform=plt.gcf().transFigure)) 

		yText = y3 + 0.015
		plt.text(.54, yText, 'AC>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.584, yText, 'AT>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.627, yText, 'CC>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.67, yText, 'CG>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.715, yText, 'CT>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.759, yText, 'GC>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.794, yText, 'TA>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.837, yText, 'TC>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.89, yText, 'TG>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.9425, yText, 'TT>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)


		if percentage:
			ytick_offest = round((y/4), 1)
			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
		else:
			# if y < 10:
			# 	if y/4 - int(y/4) > 0.5:
			# 		ytick_offest = int(y/4) + 1
			# 	else:
			# 		ytick_offest = int(y/4)
			if y < 4:
				y = 4
			#else:
			ytick_offest = int(y/4)
			if ytick_offest == 0:
				ytick_offest = 1
			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
			ylabels= [0, ytick_offest, ytick_offest*2, 
					  ytick_offest*3, ytick_offest*4]
		if y < 4:
			y = 4
		
		plt.text(0.53, 0.91+0.02, "DBS-78", fontsize=10, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
		if not percentage:
			ylabels = ['{:,}'.format(int(x)) for x in ylabels]
			if len(ylabels[-1]) > 3:
				ylabels_temp = []
				if len(ylabels[-1]) > 7:
					for label in ylabels:
						if len(label) > 7:
							ylabels_temp.append(label[0:-8] + "m")
						elif len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)

				else:
					for label in ylabels:
						if len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)
				ylabels = ylabels_temp
		# if not percentage:
		# 	ylabels = ['{:,}'.format(int(x)) for x in ylabels]
		labs = np.arange(0.44,78.44,1)
		panel7.set_xlim([0, 78])
		panel7.set_ylim([0, y])
		panel7.set_xticks(labs)
		panel7.set_yticks(ylabs)
		panel7.set_xticklabels(xlabels, rotation='vertical', fontsize=7, color='grey', fontname='Courier New', verticalalignment='top', fontweight='bold')

		panel7.set_yticklabels(ylabels, fontsize=8, color='black')
		panel7.yaxis.grid(True)
		panel7.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
		panel7.set_xlabel('')
		panel7.set_ylabel('')

		if percentage:
			panel7.set_ylabel("Percentage of Double Base Substitutions", fontsize=6, fontname="Times New Roman", weight = 'bold')
		else:
			panel7.set_ylabel("Number of Double Base Substitutions", fontsize=6, fontname="Times New Roman", weight = 'bold')
		panel7.axis('on')

		panel7.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=True,\
						   left=True, labelleft=True,\
						   right=True, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=2, colors='lightgray', width=.2)

		[i.set_color("black") for i in panel7.get_yticklabels()]
		[i.set_color("grey") for i in panel7.get_xticklabels()]


	############### plot DBS-312 ##############################################
		total_count = sum(sum(sum(tsb) for tsb in nuc.values()) for nuc in mutations_312[sample].values())
		xlabels = []

		x = 0.3
		ymax = 0
		i = 0
		colors = [[3/256,189/256,239/256], [3/256,102/256,204/256],[162/256,207/256,99/256], 
				  [1/256,102/256,1/256], [255/256,153/256,153/256], [228/256,41/256,38/256], 
				  [255/256,178/256,102/256], [255/256,128/256,1/256], [204/256,153/256,255/256], 
				  [76/256,1/256,153/256]]
		for key in mutations_312[sample]:
			muts = mutations_312[sample][key].keys()
			muts = sorted(muts)
			for seq in muts:
				xlabels.append(seq)
				if percentage:
					try:
						trans = panel8.bar(x, mutations_312[sample][key][seq][0]/total_count*100,width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
						x += 0.2
						untrans = panel8.bar(x, mutations_312[sample][key][seq][1]/total_count*100,width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
						x += .8
						if mutations_312[sample][key][seq][0]/total_count*100 > ymax:
								ymax = mutations_312[sample][key][seq][0]/total_count*100
						if mutations_312[sample][key][seq][1]/total_count*100 > ymax:
								ymax = mutations_312[sample][key][seq][1]/total_count*100
					except:
						trans = plt.bar(x, 0,width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
						untrans = plt.bar(x, 0, width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')

				else:
					trans = panel8.bar(x, mutations_312[sample][key][seq][0],width=0.2,color=[1/256,70/256,102/256],align='center', zorder=1000, label='Transcribed Strand')
					x += 0.2
					untrans = panel8.bar(x, mutations_312[sample][key][seq][1],width=0.2,color=[228/256,41/256,38/256],align='center', zorder=1000, label='Untranscribed Strand')
					x += .8
					if mutations_312[sample][key][seq][0] > ymax:
							ymax = mutations_312[sample][key][seq][0]
					if mutations_312[sample][key][seq][1] > ymax:
							ymax = mutations_312[sample][key][seq][1]
			i += 1


		y3 = .6825 + 0.02
		y = int(ymax*1.25)

		panel8.add_patch(plt.Rectangle((.785,y3), .0475, .01, facecolor=colors[0], clip_on=False, transform=plt.gcf().transFigure)) 
		panel8.add_patch(plt.Rectangle((.834,y3), .0475, .01, facecolor=colors[4], clip_on=False, transform=plt.gcf().transFigure)) 
		panel8.add_patch(plt.Rectangle((.883,y3), .0475, .01, facecolor=colors[7], clip_on=False, transform=plt.gcf().transFigure)) 
		panel8.add_patch(plt.Rectangle((.932,y3), .0475, .01, facecolor=colors[9], clip_on=False, transform=plt.gcf().transFigure)) 

		yText = y3 + .0125
		plt.text(.795, yText, 'CC>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.845, yText, 'CT>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.893, yText, 'TC>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.943, yText, 'TT>NN', fontsize=7, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		
		if y <= 4:
			y += 4

		while y%4 != 0:
			y += 1
		ytick_offest = int(y/4)

		x_shaded = 0
		panel8.add_patch(plt.Rectangle((x_shaded,0), 8.9, y, facecolor=colors[0], zorder=0, alpha = 0.25, edgecolor='grey'))
		x_shaded += 8.9
		panel8.add_patch(plt.Rectangle((x_shaded,0), 9, y, facecolor=colors[4], zorder=0, alpha = 0.25, edgecolor='grey'))
		x_shaded += 9
		panel8.add_patch(plt.Rectangle((x_shaded,0), 9, y, facecolor=colors[7], zorder=0, alpha = 0.25, edgecolor='grey'))
		x_shaded += 9
		panel8.add_patch(plt.Rectangle((x_shaded,0), 9.1, y, facecolor=colors[9], zorder=0, alpha = 0.25, edgecolor='grey'))


		if percentage:
			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
		else:

			if ytick_offest == 0:
				ytick_offest = 1
			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
			ylabels= [0, ytick_offest, ytick_offest*2, 
					  ytick_offest*3, ytick_offest*4]


		plt.text(0.78, 0.71+0.02, "DBS-186", fontsize=10, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

		if not percentage:
			ylabels = ['{:,}'.format(int(x)) for x in ylabels]
			if len(ylabels[-1]) > 3:
				ylabels_temp = []
				if len(ylabels[-1]) > 7:
					for label in ylabels:
						if len(label) > 7:
							ylabels_temp.append(label[0:-8] + "m")
						elif len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)

				else:
					for label in ylabels:
						if len(label) > 3:
							ylabels_temp.append(label[0:-4] + "k")
						else:
							ylabels_temp.append(label)
				ylabels = ylabels_temp


		labs = np.arange(0.55,36.44,1)
		panel8.set_xlim([0, 36])
		panel8.set_ylim([0, y])
		panel8.set_xticks(labs)
		panel8.set_yticks(ylabs)
		panel8.set_xticklabels(xlabels, rotation='vertical', fontsize=7, color='grey', fontname='Courier New', verticalalignment='top', fontweight='bold')

		panel8.set_yticklabels(ylabels, fontsize=8)
		panel8.yaxis.grid(True)
		panel8.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
		panel8.set_xlabel('')
		panel8.set_ylabel('')
		panel8.legend(handles=[trans, untrans], prop={'size':4})

		if percentage:
			panel8.set_ylabel("Percentage of Double Base Substitutions", fontsize=6, fontname="Times New Roman", weight = 'bold')
		else:
			panel8.set_ylabel("Number of Double Base Substitutions", fontsize=6, fontname="Times New Roman", weight = 'bold')

		panel8.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=True,\
						   left=True, labelleft=True,\
						   right=True, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=2, colors='lightgray', width=1)

		[i.set_color("black") for i in panel8.get_yticklabels()]
		[i.set_color("grey") for i in panel8.get_xticklabels()]


		panel8.set_xlim([0, 36])

		pp.savefig(plot1)
		plt.close()
	pp.close()


def main():

	# 	print("Context not supported.")
	#plotSBS("/Users/ebergstr/Desktop/Perl_tests/testCode/simulation_code_python/mutation_simulation/references/matrix/BRCA_test/BRCA_test.SBS96.all", "/Users/ebergstr/Desktop/", "test", '96', False)
	samplePortrait("/Users/ebergstr/Desktop/Mel/", "/Users/ebergstr/Desktop/Mel/output/plots/", "Mel")
	#plotDBS("/Users/ebergstr/Downloads/Biliary-AdenoCA.dinucs.csv", "/Users/ebergstr/Desktop/", "test", '78', False)

if __name__ == '__main__':
	main()
