#!/usr/bin/env python3

#This file is part of Mutational Signatures Project.

#Mutational Signatures Project: need info on project

#Copyright (C) 2018 Erik Bergstrom

#

#Mutational Signatures is free software: need to include distribtution

#rights and so forth

#

#Mutational Signatures is distributed in the hope that it will be useful,

#but WITHOUT ANY WARRANTY; without even the implied warranty of

#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

#GNU General Public License for more details [change to whatever group we should include.
 

#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu


import re
import os
import argparse
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager
from matplotlib.backends.backend_pdf import PdfPages


def plot96(matrix_path, output_path, signature, project, percentage=False):
	if 'roman' in matplotlib.font_manager.weight_dict:
		del matplotlib.font_manager.weight_dict['roman']
		matplotlib.font_manager._rebuild()

	pp = PdfPages(output_path + '96_mutations_sample_' + project + '.pdf')


	mutations = dict()
	total_count = []
	with open (matrix_path) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		for sample in samples:
			mutations[sample] = {'C>A':OrderedDict(), 'C>G':OrderedDict(), 'C>T':OrderedDict(),
								 'T>A':OrderedDict(), 'T>C':OrderedDict(), 'T>G':OrderedDict()}

		for lines in f:
			line = lines.strip().split()
			nuc = line[0]
			mut_type = line[0][2:5]
			sample_index = 1

			for sample in samples:
				if percentage:
					mutCount = float(line[sample_index])
				else:
					mutCount = int(line[sample_index])
				mutations[sample][mut_type][nuc] = mutCount
				sample_index += 1

	for sample in mutations.keys():
		total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
		plt.rcParams['axes.linewidth'] = 2
		plot1 = plt.figure(figsize=(43.93,9.92))
		plt.rc('axes', edgecolor='lightgray')
		panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
		xlabels = []
		
		x = 0.4
		ymax = 0
		colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
		i = 0
		for key in mutations[sample]:
			for seq in mutations[sample][key]:
				xlabels.append(seq[0]+seq[2]+seq[6])
				if signature:
					if percentage:
						plt.bar(x, mutations[sample][key][seq]*100,width=0.4,color=colors[i],align='center', zorder=1000)
						if mutations[sample][key][seq]*100 > ymax:
							ymax = mutations[sample][key][seq]*100
					else:	
						plt.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
						if mutations[sample][key][seq]/total_count*100 > ymax:
							ymax = mutations[sample][key][seq]/total_count*100
				else:
					plt.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq] > ymax:
							ymax = mutations[sample][key][seq]
				x += 1
			i += 1

		x = .043
		y3 = .87
		y = ymax*1.25
		y2 = y+2
		for i in range(0, 6, 1):
			panel1.add_patch(plt.Rectangle((x,y3), .15, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			x += .159

		yText = y3 + .06
		plt.text(.1, yText, 'C>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.255, yText, 'C>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.415, yText, 'C>T', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.575, yText, 'T>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.735, yText, 'T>C', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.89, yText, 'T>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

		if signature:
			ytick_offest = round((y/4), 1)
			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
		else:
			ytick_offest = int(y/4)
			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
			ylabels= [0, ytick_offest, ytick_offest*2, 
				  	  ytick_offest*3, ytick_offest*4]		

		labs = np.arange(0.375,96.375,1)

		panel1.set_xlim([0, 96])
		panel1.set_ylim([0, y])
		panel1.set_xticks(labs)
		panel1.set_yticks(ylabs)
		count = 0
		m = 0
		for i in range (0, 96, 1):
			plt.text(i/101 + .0415, .02, xlabels[i][0], fontsize=30, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			plt.text(i/101 + .0415, .044, xlabels[i][1], fontsize=30, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
			plt.text(i/101 + .0415, .071, xlabels[i][2], fontsize=30, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			count += 1
			if count == 16:
				count = 0
				m += 1	

		plt.text(0.045, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

		panel1.set_yticklabels(ylabels, fontsize=30)
		plt.gca().yaxis.grid(True)
		plt.gca().grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
		panel1.set_xlabel('')
		panel1.set_ylabel('')

		if signature:
			plt.ylabel("Mutation Percentage", fontsize=35, fontname="Times New Roman", weight = 'bold')
		else:
			plt.ylabel("Mutation Counts", fontsize=35, fontname="Times New Roman", weight = 'bold')



		panel1.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=False,\
						   left=True, labelleft=True,\
						   right=True, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors='lightgray', width=2)


		[i.set_color("black") for i in plt.gca().get_yticklabels()]

		pp.savefig(plot1)
	pp.close()


def plot192(matrix_path, output_path, signature, project, percentage=False):
	if 'roman' in matplotlib.font_manager.weight_dict:
		del matplotlib.font_manager.weight_dict['roman']
		matplotlib.font_manager._rebuild()

	pp = PdfPages(output_path + '192_mutations_sample_' + project + '.pdf')

	mutations = dict()
	with open (matrix_path) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		for sample in samples:
			mutations[sample] = {'C>A':OrderedDict(), 'C>G':OrderedDict(), 'C>T':OrderedDict(),
								 'T>A':OrderedDict(), 'T>C':OrderedDict(), 'T>G':OrderedDict()}

		for lines in f:
			line = lines.strip().split()
			nuc = line[0][2:]
			bias = line[0][0]
			if bias == 'N' or bias == 'B':
				continue
			else:
				mut_type = line[0][4:7]
				sample_index = 1

				for sample in samples:
					if percentage:
						mutCount = float(line[sample_index])
					else:
						mutCount = int(line[sample_index])
					if nuc not in mutations[sample][mut_type].keys():
						mutations[sample][mut_type][nuc] = [0,0]
					if bias == 'T':
						mutations[sample][mut_type][nuc][0] = mutCount
					else:
						mutations[sample][mut_type][nuc][1] = mutCount
					sample_index += 1

	for sample in mutations.keys():
		total_count = sum(sum(sum(tsb) for tsb in nuc.values()) for nuc in mutations[sample].values())
		plt.rcParams['axes.linewidth'] = 2
		plot1 = plt.figure(figsize=(43.93,9.92))
		plt.rc('axes', edgecolor='lightgray')
		panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
		xlabels = []
		
		x = 0.7
		ymax = 0
		colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
		i = 0
		for key in mutations[sample]:
			for seq in mutations[sample][key]:
				xlabels.append(seq[0]+seq[2]+seq[6])
				if signature:
					if percentage:
						plt.bar(x, mutations[sample][key][seq][0]/total_count*100,width=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000)
						x += 0.75
						plt.bar(x, mutations[sample][key][seq][1]/total_count*100,width=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000)
						x += .2475
						if mutations[sample][key][seq][0]/total_count*100 > ymax:
								ymax = mutations[sample][key][seq][0]/total_count*100
						elif mutations[sample][key][seq][1]/total_count*100 > ymax:
								ymax = mutations[sample][key][seq][1]/total_count*100
					else:
						plt.bar(x, mutations[sample][key][seq][0]/total_count*100,width=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000)
						x += 0.75
						plt.bar(x, mutations[sample][key][seq][1]/total_count*100,width=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000)
						x += .2475
						if mutations[sample][key][seq][0]/total_count*100 > ymax:
								ymax = mutations[sample][key][seq][0]/total_count*100
						elif mutations[sample][key][seq][1]/total_count*100 > ymax:
								ymax = mutations[sample][key][seq][1]/total_count*100

				else:
					plt.bar(x, mutations[sample][key][seq][0],width=0.75,color=[1/256,70/256,102/256],align='center', zorder=1000)
					x += 0.75
					plt.bar(x, mutations[sample][key][seq][1],width=0.75,color=[228/256,41/256,38/256],align='center', zorder=1000)
					x += .2475
					if mutations[sample][key][seq][0] > ymax:
							ymax = mutations[sample][key][seq][0]
					elif mutations[sample][key][seq][1] > ymax:
							ymax = mutations[sample][key][seq][1]
				x += 1
			i += 1

		x = .0415
		y3 = .87
		y = ymax*1.25
		x_plot = 0
		for i in range(0, 6, 1):
			panel1.add_patch(plt.Rectangle((x,y3), .155, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			panel1.add_patch(plt.Rectangle((x_plot,0), 32, ymax*1.25, facecolor=colors[i], zorder=0, alpha = 0.25, edgecolor='grey'))
			x += .1585
			x_plot += 32

		yText = y3 + .06
		plt.text(.1, yText, 'C>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.255, yText, 'C>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.415, yText, 'C>T', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.575, yText, 'T>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.735, yText, 'T>C', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.89, yText, 'T>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)


		if signature:
			ytick_offest = round((y/4),1)
			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
		else:
			ytick_offest = int(y/4)
			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
			ylabels= [0, ytick_offest, ytick_offest*2, 
				  	  ytick_offest*3, ytick_offest*4]

		labs = np.arange(0.750,192.750,1)

		panel1.set_xlim([0, 96])
		panel1.set_ylim([0, y])
		panel1.set_xticks(labs)
		panel1.set_yticks(ylabs)
		count = 0
		m = 0
		for i in range (0, 96, 1):
			plt.text(i/101 + .0415, .02, xlabels[i][0], fontsize=30, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			plt.text(i/101 + .0415, .044, xlabels[i][1], fontsize=30, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
			plt.text(i/101 + .0415, .071, xlabels[i][2], fontsize=30, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
			count += 1
			if count == 16:
				count = 0
				m += 1
		
		plt.text(0.045, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

		panel1.set_yticklabels(ylabels, fontsize=30)
		plt.gca().yaxis.grid(True)
		plt.gca().grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
		panel1.set_xlabel('')
		panel1.set_ylabel('')

		if signature:
			plt.ylabel("Mutation Percentage", fontsize=35, fontname="Times New Roman", weight = 'bold')
		else:
			plt.ylabel("Mutation Counts", fontsize=35, fontname="Times New Roman", weight = 'bold')

		panel1.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=False,\
						   left=False, labelleft=True,\
						   right=False, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors=[0.6, 0.6, 0.6])

		[i.set_color("black") for i in plt.gca().get_yticklabels()]

		pp.savefig(plot1)
	pp.close()


def plotINDEL(matrix_path, output_path, signature, project, percentage=False):
	if 'roman' in matplotlib.font_manager.weight_dict:
		del matplotlib.font_manager.weight_dict['roman']
		matplotlib.font_manager._rebuild()

	pp = PdfPages(output_path + 'INDEL_mutations_sample_' + project + '.pdf')

	indel_types = ['1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5', '1:Del:C:6'
				   '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5', '1:Del:T:6'
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

	mutations = dict()
	with open (matrix_path) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		for sample in samples:
			mutations[sample] = {'1DelC':[0,0,0,0,0,0], '1DelT':[0,0,0,0,0,0], '1InsC':[0,0,0,0,0,0], '1InsT':[0,0,0,0,0,0], 
								 '2DelR':[0,0,0,0,0,0], '3DelR':[0,0,0,0,0,0], '4DelR':[0,0,0,0,0,0], '5DelR':[0,0,0,0,0,0],
								 '2InsR':[0,0,0,0,0,0], '3InsR':[0,0,0,0,0,0], '4InsR':[0,0,0,0,0,0], '5InsR':[0,0,0,0,0,0], 
								 '2DelM':[0], '3DelM':[0,0], '4DelM':[0,0,0], '5DelM':[0,0,0,0,0]}

		for lines in f:
			line = lines.strip().split()
			categories = line[0].split(":")
			mut_type = categories[0] + categories[1] + categories[2]
			repeat_size = int(categories[3])
			if categories[2] == 'M':
				repeat_size -= 1
			sample_index = 1

			for sample in samples:
				if mut_type in mutations[sample].keys():
					if percentage:
						mutCount = float(line[sample_index])
					else:
						mutCount = int(line[sample_index])
					mutations[sample][mut_type][repeat_size] = mutCount
				else:
					continue
				sample_index += 1

	for sample in mutations.keys():
		total_count = sum(sum(nuc) for nuc in mutations[sample].values())
		plt.rcParams['axes.linewidth'] = 2
		plot1 = plt.figure(figsize=(43.93,12))
		plt.rc('axes', edgecolor='black')
		panel1 = plt.axes([0.045, 0.17, 0.92, 0.65])
		xlabels = []
		
		x = 0.4
		ymax = 0
		colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256], 
				  [253/256,202/256,181/256], [252/256,138/256,106/256], [241/256,68/256,50/256], [188/256,25/256,26/256],
				  [208/256,225/256,242/256], [148/256,196/256,223/256], [74/256,152/256,201/256], [23/256,100/256,171/256],
				  [226/256,226/256,239/256], [182/256,182/256,216/256], [134/256,131/256,189/256], [98/256,64/256,155/256]]

		i = 0
		for key in mutations[sample]:
			l = 1
			for seq in mutations[sample][key]:
				xlabels.append(l)
				if signature:
					if percentage:
						plt.bar(x, seq*100,width=0.4,color=colors[i],align='center', zorder=1000)
						if seq*100 > ymax:
							ymax = seq*100

					else:
						plt.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
						if seq/total_count*100 > ymax:
							ymax = seq/total_count*100
				else:
					plt.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq > ymax:
							ymax = seq
				x += 1
				l += 1
			i += 1

		x = .0475
		y_top = .827
		y_bottom = .114
		y = ymax*1.25
		y2 = y+2
		for i in range(0, 12, 1):
			panel1.add_patch(plt.Rectangle((x,y_top), .0595, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			panel1.add_patch(plt.Rectangle((x,y_bottom), .0595, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
			x += .0665

		panel1.add_patch(plt.Rectangle((x-.001,y_top), .006, .05, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
		panel1.add_patch(plt.Rectangle((x-.001,y_bottom), .006, .05, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
		x +=.011
		panel1.add_patch(plt.Rectangle((x,y_top), .0155, .05, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
		panel1.add_patch(plt.Rectangle((x,y_bottom), .0155, .05, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
		x += .022
		panel1.add_patch(plt.Rectangle((x,y_top), .027, .05, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
		panel1.add_patch(plt.Rectangle((x,y_bottom), .027, .05, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))		
		x += .0335
		panel1.add_patch(plt.Rectangle((x,y_top), .049, .05, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
		panel1.add_patch(plt.Rectangle((x,y_bottom), .049, .05, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))




		yText = y_top + .01
		plt.text(.072, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.1385, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
		plt.text(.205, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
		plt.text(.2715, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
		plt.text(.338, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.4045, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.471, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.5375, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		plt.text(.604, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.6705, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.737, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.8035, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		plt.text(.844, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.861, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.888, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.93, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

		yText_labels_top = yText + .075
		yText_labels_bottom = y_bottom - .03
		yText_labels_bottom_sec = yText_labels_bottom - .045

		plt.text(.08, yText_labels_top, '1bp Deletion', fontsize=40, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
		plt.text(.21, yText_labels_top, '1bp Insertion', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.375, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.64, yText_labels_top, '>1bp Insertions at Repeats\n       (Deletion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.85, yText_labels_top, 'Mircohomology\n  (Deletion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		plt.text(.058, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=35, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
		plt.text(.19, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.39, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.65, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		plt.text(.85, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x = .0477
		for i in range (0, 8, 1):
			if i != 2 and i != 3:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

			x += .0665

		for i in range (0, 4, 1):
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
			x += .0665

		plt.text(x, yText_labels_bottom, '1', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .011
		plt.text(x, yText_labels_bottom, '1  2', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .022
		plt.text(x, yText_labels_bottom, '1  2  3', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .0335
		plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		if signature:
			ytick_offest = round((y/4), 1)
			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]

		else:
			ytick_offest = int(y/4)
			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
			ylabels= [0, ytick_offest, ytick_offest*2, 
				  	  ytick_offest*3, ytick_offest*4]


		labs = np.arange(0.375,83.375,1)

		panel1.set_xlim([0, 83])
		panel1.set_ylim([0, y])
		panel1.set_xticks(labs)
		panel1.set_yticks(ylabs)	

		plt.text(0.045, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

		panel1.set_yticklabels(ylabels, fontsize=30)
		plt.gca().yaxis.grid(True)
		plt.gca().grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
		panel1.set_xlabel('')
		panel1.set_ylabel('')

		if signature:
			plt.ylabel("Mutation Percentage", fontsize=35, fontname="Times New Roman", weight = 'bold')
		else:
			plt.ylabel("Mutation Counts", fontsize=35, fontname="Times New Roman", weight = 'bold')

		panel1.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=False,\
						   left=False, labelleft=True,\
						   right=False, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors='gray', width=2)

		[i.set_color("black") for i in plt.gca().get_yticklabels()]

		pp.savefig(plot1)
	pp.close()


def plotDINUC(matrix_path, output_path, signature, project, percentage=False):
	if 'roman' in matplotlib.font_manager.weight_dict:
		del matplotlib.font_manager.weight_dict['roman']
		matplotlib.font_manager._rebuild()

	pp = PdfPages(output_path + 'DINUC_mutations_sample_' + project + '.pdf')

	dinucs = ['TT>GG','TT>CG','TT>AG','TT>GC','TT>CC','TT>AC','TT>GA','TT>CA','TT>AA','AC>CA','AC>CG','AC>CT','AC>GA',
			  'AC>GG','AC>GT','AC>TA','AC>TG','AC>TT','CT>AA','CT>AC','CT>AG','CT>GA','CT>GC','CT>GG','CT>TG','CT>TC',
			  'CT>TA','AT>CA','AT>CC','AT>CG','AT>GA','AT>GC','AT>TA','TG>GT','TG>CT','TG>AT','TG>GC','TG>CC','TG>AC',
			  'TG>GA','TG>CA','TG>AA','CC>AA','CC>AG','CC>AT','CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT','CG>AT',
			  'CG>GC','CG>GT','CG>TC','CG>TA','CG>TT','TC>GT','TC>CT','TC>AT','TC>GG','TC>CG','TC>AG','TC>GA','TC>CA',
			  'TC>AA','GC>AA','GC>AG','GC>AT','GC>CA','GC>CG','GC>TA','TA>GT','TA>CT','TA>AT','TA>GG','TA>CG','TA>GC']

	revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
	mutations = dict()
	with open (matrix_path) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		for sample in samples:
			mutations[sample] = {'AC':OrderedDict(), 'AT':OrderedDict(), 'CC':OrderedDict(),
								 'CG':OrderedDict(), 'CT':OrderedDict(), 'GC':OrderedDict(),
								 'TA':OrderedDict(), 'TC':OrderedDict(), 'TG':OrderedDict(),
								 'TT':OrderedDict()}

		for lines in f:
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
				else:
					mutCount = int(line[sample_index])
				mutations[sample][mut_type][nuc] = mutCount
				sample_index += 1

	for sample in mutations.keys():
		total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
		plt.rcParams['axes.linewidth'] = 4
		plot1 = plt.figure(figsize=(43.93,9.92))
		plt.rc('axes', edgecolor='grey')
		panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
		xlabels = []
		
		x = 0.4
		ymax = 0
		colors = [[3/256,189/256,239/256], [3/256,102/256,204/256],[162/256,207/256,99/256], 
				  [1/256,102/256,1/256], [255/256,153/256,153/256], [228/256,41/256,38/256], 
				  [255/256,178/256,102/256], [255/256,128/256,1/256], [204/256,153/256,255/256], 
				  [76/256,1/256,153/256]]
		i = 0
		for key in mutations[sample]:
			muts = mutations[sample][key].keys()
			muts = sorted(muts)
			for seq in muts:
				xlabels.append(seq)
				if signature:
					if percentage:
						plt.bar(x, mutations[sample][key][seq]*100,width=0.4,color=colors[i],align='center', zorder=1000)
						if mutations[sample][key][seq]*100 > ymax:
								ymax = mutations[sample][key][seq]*100	
					else:
						plt.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
						if mutations[sample][key][seq]/total_count*100 > ymax:
								ymax = mutations[sample][key][seq]/total_count*100					
				else:
					plt.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq] > ymax:
							ymax = mutations[sample][key][seq]
				x += 1
			i += 1

		x = .043
		y3 = .87
		y = ymax*1.25
		y2 = y+2
		i = 0
		panel1.add_patch(plt.Rectangle((.043,y3), .101, .05, facecolor=colors[0], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((.151,y3), .067, .05, facecolor=colors[1], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((.225,y3), .102, .05, facecolor=colors[2], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((.334,y3), .067, .05, facecolor=colors[3], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((.408,y3), .102, .05, facecolor=colors[4], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((.517,y3), .067, .05, facecolor=colors[5], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((.591,y3), .067, .05, facecolor=colors[6], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((.665,y3), .102, .05, facecolor=colors[7], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((.774,y3), .102, .05, facecolor=colors[8], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((.883,y3), .102, .05, facecolor=colors[9], clip_on=False, transform=plt.gcf().transFigure)) 

		yText = y3 + .06
		plt.text(.07, yText, 'AC>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.163, yText, 'AT>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.255, yText, 'CC>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.345, yText, 'CG>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.435, yText, 'CT>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.527, yText, 'GC>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.6, yText, 'TA>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.69, yText, 'TC>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.8, yText, 'TG>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
		plt.text(.915, yText, 'TT>NN', fontsize=40, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)


		if signature:
			ytick_offest = round((y/4), 1)
			ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
			ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
					  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
		else:
			if y < 10:
				if y/4 - int(y/4) > 0.5:
					ytick_offest = int(y/4) + 1
				else:
					ytick_offest = int(y/4)

			else:
				ytick_offest = int(y/4)
			if ytick_offest == 0:
				ytick_offest = 1
			ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
			ylabels= [0, ytick_offest, ytick_offest*2, 
				  	  ytick_offest*3, ytick_offest*4]
		
		
		plt.text(0.045, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)


		labs = np.arange(0.44,78.44,1)
		panel1.set_xlim([0, 78])
		panel1.set_ylim([0, y])
		panel1.set_xticks(labs)
		panel1.set_yticks(ylabs)
		panel1.set_xticklabels(xlabels, rotation='vertical', fontsize=30, color='grey', fontname='Courier New', verticalalignment='top', fontweight='bold')

		panel1.set_yticklabels(ylabels, fontsize=25)
		plt.gca().yaxis.grid(True)
		plt.gca().grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
		panel1.set_xlabel('')
		panel1.set_ylabel('')

		if signature:
			plt.ylabel("Mutation Percentage", fontsize=35, fontname="Times New Roman", weight = 'bold')
		else:
			plt.ylabel("Mutation Counts", fontsize=35, fontname="Times New Roman", weight = 'bold')

		panel1.tick_params(axis='both',which='both',\
						   bottom=False, labelbottom=True,\
						   left=True, labelleft=True,\
						   right=True, labelright=False,\
						   top=False, labeltop=False,\
						   direction='in', length=25, colors='lightgray', width=2)

		[i.set_color("black") for i in plt.gca().get_yticklabels()]
		[i.set_color("grey") for i in plt.gca().get_xticklabels()]

		pp.savefig(plot1)
	pp.close()


def plot96_single(matrix_path, sample, signature, project, percentage=False):
	if 'roman' in matplotlib.font_manager.weight_dict:
		del matplotlib.font_manager.weight_dict['roman']
		matplotlib.font_manager._rebuild()

	#pp = PdfPages(output_path + '96_mutations_sample_' + project + '.pdf')


	mutations = dict()
	total_count = []
	with open (matrix_path) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		sample_index = samples.index(sample) + 1
		mutations[sample] = {'C>A':OrderedDict(), 'C>G':OrderedDict(), 'C>T':OrderedDict(),
							 'T>A':OrderedDict(), 'T>C':OrderedDict(), 'T>G':OrderedDict()}

		for lines in f:
			line = lines.strip().split()
			nuc = line[0]
			mut_type = line[0][2:5]

			if percentage:
				mutCount = float(line[sample_index])
			else:
				mutCount = int(line[sample_index])
			mutations[sample][mut_type][nuc] = mutCount

	total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
	plt.rcParams['axes.linewidth'] = 2
	plot1 = plt.figure(figsize=(43.93,9.92))
	plt.rc('axes', edgecolor='lightgray')
	panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
	xlabels = []
	
	x = 0.4
	ymax = 0
	colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
	i = 0
	for key in mutations[sample]:
		for seq in mutations[sample][key]:
			xlabels.append(seq[0]+seq[2]+seq[6])
			if signature:
				if percentage:
					plt.bar(x, mutations[sample][key][seq]*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]*100 > ymax:
						ymax = mutations[sample][key][seq]*100
				else:	
					plt.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]/total_count*100 > ymax:
						ymax = mutations[sample][key][seq]/total_count*100
			else:
				plt.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
				if mutations[sample][key][seq] > ymax:
						ymax = mutations[sample][key][seq]
			x += 1
		i += 1

	x = .043
	y3 = .87
	y = ymax*1.25
	y2 = y+2
	for i in range(0, 6, 1):
		panel1.add_patch(plt.Rectangle((x,y3), .15, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .159

	yText = y3 + .06
	plt.text(.1, yText, 'C>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.255, yText, 'C>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.415, yText, 'C>T', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.575, yText, 'T>A', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.735, yText, 'T>C', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.89, yText, 'T>G', fontsize=55, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

	if signature:
		ytick_offest = round((y/4), 1)
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
	else:
		ytick_offest = int(y/4)
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]		

	labs = np.arange(0.375,96.375,1)

	panel1.set_xlim([0, 96])
	panel1.set_ylim([0, y])
	panel1.set_xticks(labs)
	panel1.set_yticks(ylabs)
	count = 0
	m = 0
	for i in range (0, 96, 1):
		plt.text(i/101 + .0415, .02, xlabels[i][0], fontsize=30, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		plt.text(i/101 + .0415, .044, xlabels[i][1], fontsize=30, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
		plt.text(i/101 + .0415, .071, xlabels[i][2], fontsize=30, color='gray', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		count += 1
		if count == 16:
			count = 0
			m += 1	

	plt.text(0.045, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

	panel1.set_yticklabels(ylabels, fontsize=30)
	plt.gca().yaxis.grid(True)
	plt.gca().grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
	panel1.set_xlabel('')
	panel1.set_ylabel('')

	if signature:
		plt.ylabel("Mutation Percentage", fontsize=35, fontname="Times New Roman", weight = 'bold')
	else:
		plt.ylabel("Mutation Counts", fontsize=35, fontname="Times New Roman", weight = 'bold')



	panel1.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=True, labelleft=True,\
					   right=True, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=25, colors='lightgray', width=2)


	[i.set_color("black") for i in plt.gca().get_yticklabels()]

	return(plot1)
#	pp.savefig(plot1)
#pp.close()

def main():
	signature = False
	parser = argparse.ArgumentParser(description="Provide the necessary arguments to install the reference files.")
	parser.add_argument("--matrix", "-m", help="Provide the matrix file name")
	parser.add_argument("--project", "-p", help="Provide the unique project name.")
	parser.add_argument("--context", "-c", help="Provide the unique project name.")
	parser.add_argument("-s", "--signature", help="Optional Parameter: Create the plots on a signature bases.", action='store_true')
	parser.add_argument("-sp", "--signaturePercentages", help="Optional Parameter: Create the plots on a signature bases using NMF results.", action='store_true')
	args = parser.parse_args()

	matrix = args.matrix
	project = args.project
	context = args.context
	percentage = args.signaturePercentages

	if args.signature:
		signature = True

	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)
	matrix_path = ref_dir + "/references/matrix/" + project + "/" + matrix
	output_path = ref_dir + "/plots/" + project  + "/"
	if not os.path.exists(output_path):
		os.system("mkdir " + output_path)

	if context == '96':
		plot96(matrix_path, output_path, signature, project, percentage)
	elif context == '192':
		plot192(matrix_path, output_path, signature, project, percentage)
	elif context == 'INDEL':
		plotINDEL(matrix_path, output_path, signature, project, percentage)
	elif context == 'DINUC':
		plotDINUC(matrix_path, output_path, signature, project, percentage)
	else:
		print("Context not supported.")

if __name__ == '__main__':
	main()
