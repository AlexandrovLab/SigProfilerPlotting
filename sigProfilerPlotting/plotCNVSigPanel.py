# import CNV_plotting
# import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib.font_manager
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib.patches as patches
# import pandas as pd
# import numpy as np
# from matplotlib.ticker import FixedLocator, FixedFormatter
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# import matplotlib.lines as lines
# from sklearn.metrics.pairwise import cosine_similarity
# from scipy.spatial.distance import euclidean, pdist, squareform, cosine
# import matplotlib.transforms as transforms
# import seaborn as sns
# from collections import defaultdict
# import string
# import sys
# import os
# plt.style.use('ggplot')
# plt.rcParams['axes.facecolor'] = 'white'

# def plotCNVPanel(sig_matrix):

#     df = pd.read_csv(sig_matrix, sep=None, engine='python')
#     project = sig_matrix.split("/")[-1].split('.')[0]
#     signatures = list(df)[1:]
#     n = len(signatures)
#     print ("There are " + str(n) + " signatures for " + project)

#     #AUTOMATICALLY FIGURE OUT PANEL DIMENSIONS BASED ON NUMBER OF SIGNATURES
#     if n==1:
#         num_rows=1
#         num_columns=1
#     elif n==2:
#         num_rows=1
#         num_columns=2
#     else: #n>2
#         num_rows = round((-2.78096317) * (10**-16) * (n**2) + (4.953560372*10**-1 * n) + (3.034055728*(10**-1)))  # of rows depends on number of signatures
#         num_columns = 2 #always have 2 columns

#     if n ==1:
#         fig_width=20
#         fig_height=6*n
#     else:
#         fig_width=40
#         fig_height=6*n

#     fig = plt.figure(figsize=(fig_width, fig_height)) #adjust total size of figure here
#     alpha_dict = dict(enumerate(string.ascii_lowercase))
#     for i, (col, signature) in enumerate(zip(df.columns[1:], signatures)):
#         counts = list(df[col])
#         counts = [(x/sum(counts))*100 for x in counts]
#         title = 'CNV48' + alpha_dict[i]

#         ax = fig.add_subplot(num_rows, num_columns, i+1) #ADJUST DEPENDING ON NUMBER OF SIGNATURES

#         super_class = ['Het', 'LOH', "Hom del"]
#         hom_del_class = ['0-100kb', '100kb-1Mb', '>1Mb']
#         loh_subclass = ["1", '2', '3-4', '5-8', '9+']
#         het_sub_class = ['2', '3-4', '5-8', '9+']
#         x_labels = ['0-100kb', '100kb-1Mb', '1Mb-10Mb', '10Mb-40Mb','>40Mb']
#         color_mapping = {'9+':{'>40Mb':"deeppink", '10Mb-40Mb':"hotpink", '1Mb-10Mb':"palevioletred", '100kb-1Mb':"lightpink", '0-100kb':"lavenderblush"},
#                      '5-8':{'>40Mb':"saddlebrown", '10Mb-40Mb':"sienna", '1Mb-10Mb': "peru", '100kb-1Mb':"sandybrown", '0-100kb':"linen"},
#                      '3-4':{'>40Mb': "rebeccapurple", '10Mb-40Mb':"blueviolet", '1Mb-10Mb':"mediumorchid", '100kb-1Mb':"plum", '0-100kb':"thistle"},
#                      '2':{'>40Mb':"olive", '10Mb-40Mb':"olivedrab", '1Mb-10Mb':"yellowgreen", '100kb-1Mb':"lawngreen", '0-100kb':"greenyellow"},
#                      '1':{'>40Mb':"dimgray", '10Mb-40Mb':"darkgrey", '1Mb-10Mb':"silver", '100kb-1Mb':"lightgray", '0-100kb':"whitesmoke"}}

#         hom_del_color_mapping = {'0-100kb':"darkblue", '100kb-1Mb':"mediumblue", '>1Mb':"cornflowerblue"}
#         hom_del_color_mapping = {'0-100kb':"cornflowerblue", '100kb-1Mb':"mediumblue", '>1Mb':"darkblue"}
#         patch_colors = ['green', 'purple', 'darkorange', 'fuchsia', 'slategrey', 'green', 'purple', 'darkorange', 'fuchsia', 'slateblue']

#         N=48
#         ticks = np.arange(N)
#         width = 0.27
#         xticks = []
#         i = -1 #used to distinguish first bar from the rest

#         fig, ax = plt.subplots(figsize=(16,10))

#         for count, label in zip(counts, labels):
#             categories = label.split(':')
#             cnv_class = categories[0]
#             size_class = categories[2]

#             #hom del has different color scheme and size classification
#             hom_del = False
#             if categories[1] == "homdel":
#                 hom_del = True
#             i += 1 #position of bar
#             if i == 0: #very first bar
#                 if count > 0:
#                     if hom_del:
#                         ax.bar(ticks[i], count, color=hom_del_color_mapping[size_class], edgecolor='black')
#                     else:
#                         ax.bar(ticks[i], count, color=color_mapping[cnv_class][size_class], edgecolor='black')
#                 xticks.append(ticks[i])
#             else: #all the bars besides the first bar
#                 if count > 0:
#                     if hom_del:
#                         ax.bar(ticks[i], count, color=hom_del_color_mapping[size_class], edgecolor='black')
#                     else:
#                         ax.bar(ticks[i], count, color=color_mapping[cnv_class][size_class], edgecolor='black')

#                 xticks.append(ticks[i])

#         #ADD PATCHES AND TEXT
#         patch_height = 0.2
#         patch_width = 5
#         left_edge = 0.151 #placement of left edge of patch
#         y_pos = 0.95#placement of patch on y-axis
#         text_height = 0.96

#         categories = het_sub_class + loh_subclass + ['Hom' + '\n' + 'Del']

#         trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

#         patch_locs = np.arange(0, 45, 5) #position of patches in data coordinates
#         line_locs = [] # for recording positions of top evel patches and separation lines
#         for i, loc in enumerate(patch_locs): #add 10 patches
#             ax.add_patch(plt.Rectangle((loc-0.5, 1), patch_width, patch_height, clip_on=False, facecolor=patch_colors[i], transform=trans))
#             plt.text(loc+0.5, 1.1, categories[i].capitalize(), fontsize=36, fontname='Times New Roman', fontweight='bold', color='white', transform=trans)
#             if i == 4:
#                 ax.axvline(x=loc-0.5, color='black', linewidth=2)
#                 line_locs.append(loc-0.5)

#         #add final patch for hom del
#         y_text = 0.095
#         ax.add_patch(plt.Rectangle((45-0.5, 1), 4, patch_height*2, clip_on=False, facecolor='blue', transform=trans)) #for hom del
#         ax.axvline(x=45-0.5, color='black', linewidth=2)
#         plt.text(45-0.5, 1.2+.05, "Hom-", fontsize=36, fontname='Times New Roman', fontweight='bold', color='white', transform=trans)
#         plt.text(45, 1.05, "Del", fontsize=36, fontname='Times New Roman', fontweight='bold', color='white', transform=trans)

#         #manually add top level patches(het and LOH) and text inside patches
#         ax.add_patch(plt.Rectangle((-0.5, 1.2), patch_width*4, patch_height, clip_on=False, facecolor='gray', transform=trans))
#         plt.text(8.5, 1.2+.05, "Het", fontsize=42, fontname='Times New Roman', fontweight='bold', color='white', transform=trans)
#         ax.add_patch(plt.Rectangle((line_locs[0], 1.2), patch_width*5, patch_height, clip_on=False, facecolor='black', transform=trans))
#         plt.text(line_locs[0]+10, 1.2+.05, "LOH", fontsize=42, fontname='Times New Roman', fontweight='bold', color='white', transform=trans)

#         ax.set_xticks(xticks);
#         ax.set_xticklabels(x_labels * 9 + hom_del_class, rotation=90, weight="bold", fontsize = 16);
#         ax.tick_params(labelleft=True, left=False, bottom=False)
#         ax.tick_params(axis='y', which='major', pad=0, labelsize=30)

#         #hide y-axis ticks and labels
#         # plt.gca().get_yaxis().set_ticks([])
#         # plt.gca().get_yaxis().set_ticklabels([])

#         #y-axis label
#         if aggregate:
#             ax.set_ylabel("# of events per sample", fontsize=24, fontname="Times New Roman", weight = 'bold', labelpad = 8)
#         elif percentage:
#             ax.set_ylabel("Percentage(%)", fontsize=24, fontname="Times New Roman", weight = 'bold', labelpad = 8)
#             ax.yaxis.labelpad = 1
#         else:
#             ax.set_ylabel("# of events", fontsize=24, fontname="Times New Roman", weight = 'bold', labelpad = 8)

#         #TITLE
#         if not aggregate:
#             plt.text(0, 0.90, sample, fontsize=20, fontname='Times New Roman', fontweight='bold', color='black', transform=trans)
#         else:
#             plt.text(0, 0.90, project, fontsize=20, fontname='Times New Roman', fontweight='bold', color='black', transform=trans)

#         plt.tight_layout()
#         plt.rcParams["font.weight"] = "bold"

#     #SAVE FILE
#     fig.tight_layout()
#     if n < 10:
#         fig.savefig(output_path + 'CNV48_S' + str(n) + '_signature_panel.png', dpi=600, bbox_inches='tight')
#         print("Saved plot to " + output_path + 'CNV48_S' + str(n) + '_signature_panel.png')
#     else:
#         fig.savefig(output_path + 'CNV48_S' + str(n) + '_signature_panel.png', dpi=300, bbox_inches='tight')
#         print("Saved plot to " + output_path + 'CNV48_S' + str(n) + '_signature_panel.png')
