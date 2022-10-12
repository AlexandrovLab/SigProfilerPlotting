# #!/usr/bin/env python3

# #Author: Erik Bergstrom

# #Contact: ebergstr@eng.ucsd.edu

# from bdb import set_trace
# from copy import deepcopy
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib.font_manager
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib.patches as mplpatches
# import matplotlib.ticker as ticker
# from matplotlib.ticker import LinearLocator
# import matplotlib.lines as lines
# import matplotlib.transforms as transforms
# import re
# import os
# import sys
# import argparse
# from collections import OrderedDict
# import pandas as pd
# import numpy as np
# import io
# import string
# import warnings
# import pickle
# import sigProfilerPlotting as spplt
# import pdb
# import itertools,time
# import sklearn
# from sklearn.preprocessing import LabelEncoder


# def getIDtemp(output_path,project):
#     # pp = PdfPages(output_path + 'ID83_template_' + project + '.pdf')
#     plt.rcParams['axes.linewidth'] = 2
#     plot1 = plt.figure(figsize=(43.93,12))
#     plt.rc('axes', edgecolor='black')
#     panel1 = plt.axes([0.045, 0.17, 0.92, 0.65])
#     xlabels = []

#     x = 0.4
#     ymax = 0
#     colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256],
#                 [253/256,202/256,181/256], [252/256,138/256,106/256], [241/256,68/256,50/256], [188/256,25/256,26/256],
#                 [208/256,225/256,242/256], [148/256,196/256,223/256], [74/256,152/256,201/256], [23/256,100/256,171/256],
#                 [226/256,226/256,239/256], [182/256,182/256,216/256], [134/256,131/256,189/256], [98/256,64/256,155/256]]

#     x = .0475
#     y_top = .827
#     y_bottom = .114
#     y = int(ymax*1.25)
#     y2 = y+2
#     for i in range(0, 12, 1):
#         panel1.add_patch(plt.Rectangle((x,y_top), .0595, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
#         panel1.add_patch(plt.Rectangle((x,y_bottom), .0595, .05, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure))
#         x += .0665

#     panel1.add_patch(plt.Rectangle((x-.001,y_top), .006, .05, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
#     panel1.add_patch(plt.Rectangle((x-.001,y_bottom), .006, .05, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
#     x +=.011
#     panel1.add_patch(plt.Rectangle((x,y_top), .0155, .05, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
#     panel1.add_patch(plt.Rectangle((x,y_bottom), .0155, .05, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
#     x += .022
#     panel1.add_patch(plt.Rectangle((x,y_top), .027, .05, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
#     panel1.add_patch(plt.Rectangle((x,y_bottom), .027, .05, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
#     x += .0335
#     panel1.add_patch(plt.Rectangle((x,y_top), .049, .05, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
#     panel1.add_patch(plt.Rectangle((x,y_bottom), .049, .05, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))




#     yText = y_top + .01
#     plt.text(.072, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
#     plt.text(.1385, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
#     plt.text(.205, yText, 'C', fontsize=40, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
#     plt.text(.2715, yText, 'T', fontsize=40, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
#     plt.text(.338, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.4045, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.471, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.5375, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
#     plt.text(.604, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.6705, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.737, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.8035, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
#     plt.text(.844, yText, '2', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.861, yText, '3', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.888, yText, '4', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.93, yText, '5+', fontsize=40, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

#     yText_labels_top = yText + .075
#     yText_labels_bottom = y_bottom - .03
#     yText_labels_bottom_sec = yText_labels_bottom - .045

#     plt.text(.08, yText_labels_top, '1bp Deletion', fontsize=40, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
#     plt.text(.21, yText_labels_top, '1bp Insertion', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.375, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.64, yText_labels_top, '>1bp Insertion at Repeats\n       (Insertion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.85, yText_labels_top, ' Microhomology\n(Deletion Length)', fontsize=40, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

#     plt.text(.058, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=35, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
#     plt.text(.19, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.39, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.65, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     plt.text(.85, yText_labels_bottom_sec, 'Microhomology Length', fontsize=35, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

#     x = .0477
#     for i in range (0, 8, 1):
#         if i != 2 and i != 3:
#             plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#         else:
#             plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

#         x += .0665

#     for i in range (0, 4, 1):
#         plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#         x += .0665

#     plt.text(x, yText_labels_bottom, '1', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     x += .011
#     plt.text(x, yText_labels_bottom, '1  2', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     x += .022
#     plt.text(x, yText_labels_bottom, '1  2  3', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
#     x += .0335
#     plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=32, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)


#     labs = np.arange(0.375,83.375,1)

#     panel1.set_xlim([0, 83])
#     panel1.set_ylim([0, y])
#     panel1.set_xticks(labs)



#     # panel1.set_yticklabels(ylabels, fontsize=30)
#     plt.gca().yaxis.grid(True)
#     plt.gca().grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
#     panel1.set_xlabel('')
#     panel1.set_ylabel('')

#     panel1.tick_params(axis='both',which='both',\
#                         bottom=False, labelbottom=False,\
#                         left=False, labelleft=True,\
#                         right=False, labelright=False,\
#                         top=False, labeltop=False,\
#                         direction='in', length=25, colors='gray', width=2)

#     [i.set_color("black") for i in plt.gca().get_yticklabels()]
#     package_path = spplt.__path__[0]
#     path_1 =os.path.join(package_path,'templates/')
#     filename= os.path.join(path_1,'ID'+'.pkl')
#     pickle.dump(plot1, open(filename, 'wb'))
#     # pp.savefig(plot1)
#     # pp.close()



# def plotIDmod(matrix_path, output_path, project, plot_type, percentage=False, custom_text_upper=None, custom_text_middle=None, custom_text_bottom=None,savefig_format = "pdf"):
#     # if 'roman' in matplotlib.font_manager.weight_dict:
#     #   del matplotlib.font_manager.weight_dict['roman']
#     #   matplotlib.font_manager._rebuild()
#     # import pdb;pdb.set_trace()
#     plot_custom_text = False
#     sig_probs = False
#     pcawg = False
#     if plot_type == '94' or plot_type == 'ID94' or plot_type == '94ID' or plot_type == '83':
        
#         # getIDtemp(output_path,project)
#         # import pdb;pdb.set_trace()
#         if not isinstance(matrix_path, pd.DataFrame):
#             data=pd.read_csv(matrix_path,sep='\t',index_col=0)
#             data=data.dropna(axis=1, how='all')

#         if data.isnull().values.any():
#             raise ValueError("Input data contains Nans.")


#         try:


#             sample_count = 0
#             buf= io.BytesIO()
#             fig_orig=pickle.load(open(spplt.__path__[0]+'/templates/ID.pkl','rb'))
#             pickle.dump(fig_orig, buf)

#             figs={}

#             colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256],
#                             [253/256,202/256,181/256], [252/256,138/256,106/256], [241/256,68/256,50/256], [188/256,25/256,26/256],
#                             [208/256,225/256,242/256], [148/256,196/256,223/256], [74/256,152/256,201/256], [23/256,100/256,171/256],
#                             [226/256,226/256,239/256], [182/256,182/256,216/256], [134/256,131/256,189/256], [98/256,64/256,155/256]]
#             ctx = data.index
#             xlabels= [i.split(":")[0]+i.split(":")[1]+i.split(":")[2] for i in ctx.to_list()]
#             xlables_set= ['1DelC', '1DelT', '1InsC', '1InsT', '2DelR', '3DelR', '4DelR', '5DelR', '2InsR','3InsR','4InsR','5InsR', '2DelM', '3DelM','4DelM','5DelM']
#             colors_idx = deepcopy(xlabels)
#             for ii in range(0,len(xlables_set)):
#                 colors_idx=[ii if x==xlables_set[ii] else x for x in colors_idx]
            
#             # pdb.set_trace()
#             colors_flat_list = [colors[i] for i in colors_idx]

#             for sample in data.columns: #mutations.keys():
                
#                 buf.seek(0)
#                 figs[sample]=pickle.load(buf)
#                 panel1= figs[sample].axes[0]
#                 muts= data[sample].values
#                 total_count = np.sum(muts) 
#                 x = 0.4

#                 if percentage:
#                     if total_count > 0:
#                         plt.bar(np.arange(len(ctx))+x, muts/total_count*100,width=0.4,color=colors_flat_list,align='center', zorder=1000)
#                         ymax = np.max(muts/total_count*100)
#                 else:
#                     plt.bar(np.arange(len(ctx))+x, muts,width=0.4,color=colors_flat_list,align='center', zorder=1000)
#                     ymax = np.max(muts)

#                 x = .0475
#                 y_top = .827
#                 y_bottom = .114
#                 y = int(ymax*1.25)
#                 y2 = y+2

#                 if y <= 4:
#                     y += 4

#                 while y%4 != 0:
#                     y += 1
#                 ytick_offest = int(y/4)

#                 if percentage:
#                     ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
#                     ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%",
#                                 str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
#                 else:
#                     ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
#                     ylabels= [0, ytick_offest, ytick_offest*2,
#                                 ytick_offest*3, ytick_offest*4]

#                 labs = np.arange(0.375,83.375,1)

#                 if not percentage:
#                     ylabels=spplt.getylabels(ylabels)

#                 panel1.set_xlim([0, 83])
#                 panel1.set_ylim([0, y])
#                 panel1.set_xticks(labs)
#                 panel1.set_yticks(ylabs)
#                 # import pdb;pdb.set_trace()
#                 if sig_probs:
#                     plt.text(0.0475, 0.75, sample, fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)
#                 else:
#                     plt.text(0.0475, 0.75, sample + ": " + "{:,}".format(int(total_count)) + " indels", fontsize=60, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure)

#                 custom_text_upper_plot = ''
#                 try:
#                     custom_text_upper[sample_count]
#                 except:
#                     custom_text_upper = False
#                 try:
#                     custom_text_middle[sample_count]
#                 except:
#                     custom_text_middle = False
#                 try:
#                     custom_text_bottom[sample_count]
#                 except:
#                     custom_text_bottom = False

#                 if custom_text_upper:
#                     plot_custom_text = True
#                     if len(custom_text_upper[sample_count]) > 40:
#                         print("To add a custom text, please limit the string to <40 characters including spaces.")
#                         plot_custom_text = False
#                 if custom_text_middle:
#                     if len(custom_text_middle[sample_count]) > 40:
#                         print("To add a custom text, please limit the string to <40 characters including spaces.")
#                         plot_custom_text = False

#                 if plot_custom_text:
#                     x_pos_custom = 0.95
#                     if custom_text_upper and custom_text_middle:
#                         custom_text_upper_plot = custom_text_upper[sample_count] + "\n" + custom_text_middle[sample_count]
#                         if custom_text_bottom:
#                             custom_text_upper_plot += "\n" + custom_text_bottom[sample_count]

#                     if custom_text_upper and not custom_text_middle:
#                         custom_text_upper_plot = custom_text_upper[sample_count]
#                         panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

#                     elif custom_text_upper and custom_text_middle:
#                         if not custom_text_bottom:
#                             panel1.text(x_pos_custom, 0.72, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')
#                         else:
#                             panel1.text(x_pos_custom, 0.68, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')

#                     elif not custom_text_upper and custom_text_middle:
#                         custom_text_upper_plot = custom_text_middle[sample_count]
#                         panel1.text(x_pos_custom, 0.78, custom_text_upper_plot, fontsize=35, weight='bold', color='black', fontname= "Arial", transform=plt.gcf().transFigure, ha='right')




#                 panel1.set_yticklabels(ylabels, fontsize=30)
#                 plt.gca().yaxis.grid(True)
#                 plt.gca().grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
#                 panel1.set_xlabel('')
#                 panel1.set_ylabel('')

#                 if percentage:
#                     plt.ylabel("Percentage of Indels", fontsize=35, fontname="Times New Roman", weight = 'bold')
#                 else:
#                     plt.ylabel("Number of Indels", fontsize=35, fontname="Times New Roman", weight = 'bold')

#                 panel1.tick_params(axis='both',which='both',\
#                                     bottom=False, labelbottom=False,\
#                                     left=False, labelleft=True,\
#                                     right=False, labelright=False,\
#                                     top=False, labeltop=False,\
#                                     direction='in', length=25, colors='gray', width=2)

#                 [i.set_color("black") for i in plt.gca().get_yticklabels()]
#                 # import pdb;pdb.set_trace()
#                 # pp.savefig(plot1)
#                 # plt.close()
#                 sample_count += 1
#             # pp.close()
#             if savefig_format == "pdf":
#                 pp = PdfPages(output_path + 'ID_83_plots_' + project + '.pdf') # PdfPages(output_path + 'SBS_96_plots_' + project + '.pdf')
        
            
#             if savefig_format == "pdf":
#                 for fig in figs:
#                     figs[fig].savefig(pp, format='pdf')
#                 pp.close()
                
#             if savefig_format == "png":
#                 for fig in figs:
#                     figs[fig].savefig(output_path + 'SBS_96_plots_'+fig+'.png',dpi=100)
#             if savefig_format == "buffer_stream":
#                 buff_list={}
#                 for fig in figs:
#                     buffer2=io.BytesIO()
#                     figs[fig].savefig(buffer2,format='png')
#                     buff_list[fig]=buffer2
#                     return buff_list
#         except:
#             print("There may be an issue with the formatting of your matrix file.")
#             os.remove(output_path + 'ID_83_plots_' + project + '.pdf')
