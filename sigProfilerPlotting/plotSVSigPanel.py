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


# 560_Breast_CNV.De_Novo_Solution_Signatures_SBS48.tsv
# 560_Breast_SV.De_Novo_Solution_Signatures_SBS32.tsv
# Mutographs_ESCC_CNV.De_Novo_Solution_Signatures_SBS48.tsv
# Mutographs_ESCC_SV_De_Novo_Solution_Signatures_SBS32.tsv

#sig_matrix = '/Users/azhark/Documents/Alexandrov_Lab/matrices/Mutographs_Rearrangements_De_Novo_Solution_Signatures_SBS32.tsv' #path to signature matrix
#sig_matrix = '/Users/azhark/Documents/Alexandrov_Lab/for_Ludmil/Mutographs_ESCC_SV_De_Novo_Solution_Signatures_SBS32.tsv'
#matrices = [] #list of signature matrices

#input_path = '/Users/azhark/Documents/Alexandrov_Lab/SigProfilerExtractor/560_breast_NoNorm_100rep_0.05pen_0.9_15/SBS32/All_solutions/'
#input_path = '/Users/azhark/Documents/Alexandrov_Lab/SigProfilerExtractor/Mutographs_ESCC_NoNorm_100rep_0.05pen_0.9_15/SBS32/All_solutions/'
input_path = "/Users/khandekara2/iCloud/SigProfilerExtractor/Mutographs_ESCC_NoNorm_100rep_0.05pen_0.9_12/SBS32/All_Solution_Plots/"

#output_path = input_path + "All_Signature_Plots/"
# output_path = "/Users/khandekara2/iCloud/SigProfilerExtractor/Mutographs_ESCC_NoNorm_100rep_0.05pen_0.9_12/SBS32/All_Solution_Plots/"
#
# if not os.path.exists(output_path):
#     os.makedirs(output_path)
#     print("Output Directory '%s' created" %output_path)
# else:
#     print("Output Directory '%s' already exists so plots will be outputted there" %output_path)
#
# matrices = []
# os.chdir(input_path)
# for file in os.listdir("."):
#     if file.endswith(".txt"):
#         matrices.append(file)

#matrices = ["Sherlock-Lung_Rearrangement_Signatures_SBS32.tsv"]
# for sig_matrix in matrices:
#     print (sig_matrix)
#     df = pd.read_csv(sig_matrix, sep=None, engine='python')
#     project = sig_matrix.split("/")[-1].split('.')[0]
#     signatures = list(df)[1:]
#     n = len(signatures)
#     print ("There are " + str(n) + " signatures for " + project)

    # print (sys.argv)
    # try:
    #     fig_width = int(sys.argv[1])
    #     fig_height = int(sys.argv[2])
    #     num_rows = int(sys.argv[3])
    #     num_columns = int(sys.argv[4])
    # except:
    #     print ("Please provide the height and width of the figure and number of rows and columns of subplots")

    #AUTOMATICALLY FIGURE OUT PANEL DIMENSIONS BASED ON NUMBER OF SIGNATURES
#     if n==1:
#         num_rows=1
#         num_columns=1
#     elif n==2:
#         num_rows=1
#         num_columns=2
#     else: #n>2
#         num_rows = round((-2.78096317) * (10**-16) * (n**2) + (4.953560372*10**-1 * n) + (3.034055728*(10**-1)))  # of rows depends on number of signatures
#         num_columns = 2 #always have 2 columns
#
#     if n ==1:
#         fig_width=20
#         fig_height=6*n
#     else:
#         fig_width=40
#         fig_height=6*n
#
#     fig = plt.figure(figsize=(fig_width, fig_height)) #adjust total size of figure here
#
#     color_mapping = {'del':{'>10Mb':"deeppink", '1Mb-10Mb':"hotpink", '10-100Kb':"palevioletred", '100Kb-1Mb':"lightpink", '1-10Kb':"lavenderblush"},
#                      'tds':{'>10Mb':"saddlebrown", '1Mb-10Mb':"sienna", '10-100Kb':"peru", '100Kb-1Mb':"sandybrown", '1-10Kb':"linen"},
#                      'inv':{'>10Mb':"rebeccapurple", '1Mb-10Mb':"blueviolet", '10-100Kb':"mediumorchid", '100Kb-1Mb':"plum", '1-10Kb':"thistle"}}
#     alpha_dict = dict(enumerate(string.ascii_lowercase))
#
#     labels = df[df.columns[0]]
#     samples = list(df)[1:]
#     x_labels = ['1-10kb', '10-100kb', '100kb-1Mb', '1Mb-10Mb','>10Mb']
#     super_class = ['clustered', 'non-clustered']
#     sub_class = ['del', 'tds', 'inv', 'trans']
#     #print (num_rows)
#     for i, (col, sample) in enumerate(zip(df.columns[1:], samples)):
#
#         title = 'RS32' + sample[-1]
#         ax = fig.add_subplot(num_rows, num_columns, i+1) #ADJUST DEPENDING ON NUMBER OF SIGNATURES
#         counts = list(df[col])
#         counts = [(x/sum(counts))*100 for x in counts] #PERCENTAGE
#         #print (len(counts))
#         assert(len(counts)) == 32
#         assert(len(labels)) == 32
#         N=32
#         ticks = [0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34]
#         #print (len(ticks))
#         width = 0.27
#         xticks = []
#         i = -1 #used to distinguish first bar from the rest
#
#         #now loop through each count in a signature column
#         for count, label in zip(counts, labels):
#             categories = label.split('_')
#             if len(categories) > 2:
#                 rearrangement_class = categories[1]
#                 size_class = categories[2]
#             i += 1 #position of bar
#             #print (categories)
#             if i == 0: #very first bar
#                 if count >= 0:
#                     if len(categories) == 2: #clustered translocation or non-clustered translocation
#                         ax.bar(ticks[i], count, color="dimgray", edgecolor='black'); #translocation only has one color
#                     else:
#                         ax.bar(ticks[i], count, color=color_mapping[rearrangement_class][size_class], edgecolor='black');
#             else: #all the bars besides the first bar
#                 if count >= 0:
#                     if len(categories) == 2: #clustered translocation or non-clustered translocation
#                         ax.bar(ticks[i], count, color="dimgray", edgecolor='black'); #translocation only has one color
#                     else:
#                         ax.bar(ticks[i], count, color=color_mapping[rearrangement_class][size_class], edgecolor='black');
#             xticks.append(ticks[i])
#
#         #ADD PATCHES AND TEXT
#         patch_height = 0.2
#         patch_width = 5
#         left_edge = 0.151 #placement of left edge of patch
#         y_pos = 0.95 #placement of patch on y-axis
#         text_height = 0.96
#         patch_colors = ['maroon', 'darkorange', 'slateblue', 'green', 'maroon', 'darkorange', 'slateblue', 'green']
#         classes = sub_class + sub_class
#
#         trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
#         patch_locs = np.arange(0, 40, 5)
#         patch_locs1 = [0, 5, 10]
#         #print (patch_locs)#position of patches in data coordinates
#         line_locs = [] # for recording positions of top level patches and separation lines
#
#         for i, loc in enumerate(patch_locs1): #add first 3 patches
#             ax.add_patch(plt.Rectangle((loc-0.5, 1), patch_width, patch_height, clip_on=False, facecolor=patch_colors[i], transform=trans));
#             plt.text(loc+1, 1.05, classes[i].capitalize(), fontsize=36, fontname='Times New Roman', fontweight='bold', color='white', transform=trans);
#
#         t="trans"
#         ax.add_patch(plt.Rectangle((15-0.5, 1), 3, patch_height, clip_on=False, facecolor='dimgray', transform=trans)); #patch for trans
#         plt.text(15-0.3, 1.05, t.capitalize(), fontsize=30, fontname='Times New Roman', fontweight='bold', color='white', transform=trans);
#
#         ax.axvline(x=18-0.5, color='black', linewidth=2);
#         line_locs.append(loc-0.5)
#
#         patch_locs2 = [18, 23, 28]
#         for i, loc in enumerate(patch_locs2): #add next 3 patches
#             ax.add_patch(plt.Rectangle((loc-0.5, 1), patch_width, patch_height, clip_on=False, facecolor=patch_colors[i], transform=trans));
#             plt.text(loc+1, 1.05, classes[i].capitalize(), fontsize=36, fontname='Times New Roman', fontweight='bold', color='white', transform=trans);
#
#         ax.add_patch(plt.Rectangle((33-0.5, 1), 3, patch_height, clip_on=False, facecolor='dimgray', transform=trans)); #patch for trans
#         plt.text(33-0.3, 1.05, t.capitalize(), fontsize=30, fontname='Times New Roman', fontweight='bold', color='white', transform=trans);
#
#         #manually add top level patches(het and LOH) and text inside patches
#         ax.add_patch(plt.Rectangle((-0.5, 1.2), patch_width*4, patch_height, clip_on=False, facecolor='gray', transform=trans));
#         plt.text(6.25, 1.2+.05, "Clustered", fontsize=42, fontname='Times New Roman', fontweight='bold', color='white', transform=trans);
#         ax.add_patch(plt.Rectangle((18-0.5, 1.2), (patch_width*3)+3, patch_height, clip_on=False, facecolor='black', transform=trans));
#         plt.text(23, 1.2+.05, "Non-Clustered", fontsize=42, fontname='Times New Roman', fontweight='bold', color='white', transform=trans);
#
#         ax.set_xticks(xticks);
#         ax.set_xticklabels(x_labels * 3 + [' '] + x_labels * 3 + [' '], rotation=90, weight="bold", fontsize = 16);
#         ax.tick_params(labelleft=True, left=False, bottom=False)
#         ax.tick_params(axis='y', which='major', pad=0, labelsize=30)
#
#     #   TITLE
#         plt.text(0, 0.84, title.upper(), fontsize=30, fontname='Times New Roman', fontweight='bold', color='black', transform=trans) #title
#
#         #y-axis label
#         ax.set_ylabel("Percentage(%)", fontsize=30, fontname="Times New Roman", weight ='bold')
#         ax.yaxis.labelpad = 1
#         plt.tight_layout()
#         plt.rcParams["font.weight"] = "bold"
#
#     fig.tight_layout()
#     if n > 10: #lower resolution figure
#         fig.savefig(output_path + 'RS32_S' + str(n) + '_signature_panel.png', dpi=300, bbox_inches='tight')
#     else:
#         fig.savefig(output_path + 'RS32_S' + str(n) + '_signature_panel.png', dpi=600, bbox_inches='tight')
#
# print("Saved panel figure to " + output_path)
