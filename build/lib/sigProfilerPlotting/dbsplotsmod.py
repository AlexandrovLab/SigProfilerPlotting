from bdb import set_trace
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mplpatches
import matplotlib.ticker as ticker
from matplotlib.ticker import LinearLocator
import matplotlib.lines as lines
import matplotlib.transforms as transforms
import re
import os
import sys
import argparse
from collections import OrderedDict
import pandas as pd
import numpy as np
import io
import string
import warnings
import pickle
import sigProfilerPlotting as spplt
import pdb
import itertools,time
from sklearn.preprocessing import LabelEncoder

warnings.filterwarnings("ignore")
def getdbstemp(context='DBS78',path='DBS78.pkl'):
    # pp = PdfPages(output_path + 'DBS_78_plots_template' + '.pdf')
    plot_custom_text = False
    pcawg = False
    sig_probs = False
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

    x = .043
    y3 = .87
    y = int(ymax*1.25)
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

    if y <= 4:
        y += 4

    while y%4 != 0:
        y += 1
    ytick_offest = int(y/4)

    labs = np.arange(0.44,78.44,1)
    panel1.set_xlim([0, 78])
    panel1.set_ylim([0, y])
    panel1.set_xticks(labs)
    # panel1.set_yticks(ylabs)
    panel1.set_xticklabels(xlabels, rotation='vertical', fontsize=30, color='grey', fontname='Courier New', verticalalignment='top', fontweight='bold')

    # panel1.set_yticklabels(ylabels, fontsize=25)
    plt.gca().yaxis.grid(True)
    plt.gca().grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
    panel1.set_xlabel('')
    panel1.set_ylabel('')

    # if percentage:
    #     plt.ylabel("Percentage of Double Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')
    # else:
    #     plt.ylabel("Number of Double Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')

    panel1.tick_params(axis='both',which='both',\
                        bottom=False, labelbottom=True,\
                        left=True, labelleft=True,\
                        right=True, labelright=False,\
                        top=False, labeltop=False,\
                        direction='in', length=25, colors='lightgray', width=2)

    [i.set_color("black") for i in plt.gca().get_yticklabels()]
    [i.set_color("grey") for i in plt.gca().get_xticklabels()]
    pickle.dump(plot1, open(path, 'wb'))

    # pp.savefig(plot1)
    # plt.close()
    # pp.close()

def plotDBSmod(matrix_path, output_path, project, plot_type, percentage=False, custom_text_upper=None, custom_text_middle=None, custom_text_bottom=None):

    # context ='DBS78'
    # package_path = spplt.__path__[0]
    # install_path =os.path.join(package_path,'templates/')
    # if not os.path.exists(install_path):
    #     os.mkdir(install_path)

    # filename= os.path.join(install_path,context+'.pkl')
    # getdbstemp(context,path=filename)
    plot_custom_text = False
    pcawg = False
    sig_probs = False
    if plot_type == '78' or plot_type == '78DBS' or plot_type == 'DBS78':
        if not isinstance(matrix_path, pd.DataFrame):
            data=pd.read_csv(matrix_path,sep='\t',index_col=0)
            data=data.dropna(axis=1, how='all')

        if data.isnull().values.any():
            raise ValueError("Input data contains Nans.")
        # import pdb;
        # pdb.set_trace()
        # with open(matrix_path) as f:
        #     next(f)
        #     first_line = f.readline()
        #     first_line = first_line.strip().split()
        #     mutation_type = first_line[0]
        #     if first_line[0][2] != ">":
        #         pcawg = True
        #     if len(mutation_type) != 5 and first_line[0][2] == ">":
        #         sys.exit("The matrix does not match the correct DBS96 format. Please check you formatting and rerun this plotting function.")
        # pp = PdfPages(output_path + 'DBS_78_plots_' + project + '.pdf')

        dinucs = ['TT>GG','TT>CG','TT>AG','TT>GC','TT>CC','TT>AC','TT>GA','TT>CA','TT>AA','AC>CA','AC>CG','AC>CT','AC>GA',
                    'AC>GG','AC>GT','AC>TA','AC>TG','AC>TT','CT>AA','CT>AC','CT>AG','CT>GA','CT>GC','CT>GG','CT>TG','CT>TC',
                    'CT>TA','AT>CA','AT>CC','AT>CG','AT>GA','AT>GC','AT>TA','TG>GT','TG>CT','TG>AT','TG>GC','TG>CC','TG>AC',
                    'TG>GA','TG>CA','TG>AA','CC>AA','CC>AG','CC>AT','CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT','CG>AT',
                    'CG>GC','CG>GT','CG>TC','CG>TA','CG>TT','TC>GT','TC>CT','TC>AT','TC>GG','TC>CG','TC>AG','TC>GA','TC>CA',
                    'TC>AA','GC>AA','GC>AG','GC>AT','GC>CA','GC>CG','GC>TA','TA>GT','TA>CT','TA>AT','TA>GG','TA>CG','TA>GC']

        revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
        mutations = OrderedDict()

        
        try:
            vals=set(list(data.index))-set(dinucs)
            vals=sorted(list(vals))
            for ech in vals:
                ech_mod = ech.split('>')[0]+'>'+revcompl(ech.split('>')[1])
                data.rename(index={ech:ech_mod}, inplace=True)

            data = data.sort_index()
            ctx = data.index
            xlabels=[dn.split('>')[1] for dn in ctx]
            colors = [[3/256,189/256,239/256], [3/256,102/256,204/256],[162/256,207/256,99/256],
                    [1/256,102/256,1/256], [255/256,153/256,153/256], [228/256,41/256,38/256],
                    [255/256,178/256,102/256], [255/256,128/256,1/256], [204/256,153/256,255/256],
                    [76/256,1/256,153/256]]
            mainlist=[dn.split('>')[0] for dn in ctx]
            le = LabelEncoder()
            colors_idxs = le.fit_transform(mainlist)
            colors_flat_list = [colors[i] for i in colors_idxs]
            sample_count = 0


            buf= io.BytesIO()
            fig_orig=pickle.load(open(spplt.__path__[0]+'/templates/DBS78.pkl','rb'))
            pickle.dump(fig_orig, buf)
            figs={}


            for sample in data.columns:
                buf.seek(0)
                figs[sample]=pickle.load(buf)
                panel1= figs[sample].axes[0]    
                total_count =  np.sum(data[sample].values)#sum(sum(nuc.values()) for nuc in mutations[sample].values())

                x = 0.4
                muts = data[sample].values
                if percentage:
                    if total_count > 0:
                        plt.bar(np.asarray(range(len(ctx)))+x,muts/total_count*100,width=0.4,color=colors_flat_list,align='center', zorder=1000)
                        ymax = np.max(muts/total_count*100)
                else:

                    plt.bar(np.asarray(range(len(ctx)))+x,muts,width=0.4,color=colors_flat_list,align='center', zorder=1000)
                    ymax = np.max(muts)
                # for i in range(len(xlabels)):
                #     print(xlabels[i],muts[i])

                x = .043
                y3 = .87
                y = int(ymax*1.25)
                y2 = y+2
                i = 0
 
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


                if sig_probs:
                    plt.text(0.045, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
                else:
                    plt.text(0.045, 0.75, sample + ": " + "{:,}".format(int(total_count)) + " double subs", fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

                custom_text_upper_plot = ''
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_middle[sample_count]
                except:
                    custom_text_middle = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print("To add a custom text, please limit the string to <40 characters including spaces.")
                        plot_custom_text = False
                if custom_text_bottom:
                    if len(custom_text_bottom[sample_count]) > 40:
                        print("To add a custom text, please limit the string to <40 characters including spaces.")
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.98
                    if custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = custom_text_upper[sample_count] + "\n" + custom_text_middle[sample_count]
                        if custom_text_bottom:
                            custom_text_upper_plot += "\n" + custom_text_bottom[sample_count]

                    if custom_text_upper and not custom_text_middle:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

                    elif custom_text_upper and custom_text_middle:
                        if not custom_text_bottom:
                            panel1.text(x_pos_custom, 0.75, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')
                        else:
                            panel1.text(x_pos_custom, 0.7, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

                    elif not custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = custom_text_middle[sample_count]
                        panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')


                if not percentage:
                    ylabels =spplt.getylabels(ylabels)

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

                if percentage:
                    plt.ylabel("Percentage of Double Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')
                else:
                    plt.ylabel("Number of Double Base Substitutions", fontsize=35, fontname="Times New Roman", weight = 'bold')

                panel1.tick_params(axis='both',which='both',\
                                    bottom=False, labelbottom=True,\
                                    left=True, labelleft=True,\
                                    right=True, labelright=False,\
                                    top=False, labeltop=False,\
                                    direction='in', length=25, colors='lightgray', width=2)

                [i.set_color("black") for i in plt.gca().get_yticklabels()]
                [i.set_color("grey") for i in plt.gca().get_xticklabels()]

                # pp.savefig(plot1)
                # plt.close()
                # sample_count += 1
            pp = PdfPages(output_path + 'DBS_78_plots_' + project + '.pdf')
            for fig in figs:
                figs[fig].savefig(pp, format='pdf')
            pp.close()
            
        except:
            print("There may be an issue with the formatting of your matrix file.")
            os.remove(output_path + 'DBS_78_plots_' + project + '.pdf')