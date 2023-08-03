import os

import sigProfilerPlotting as sigPlt

matrix_path = "input/examples/samples/"
output_path = "output/examples/"

if not os.path.exists(output_path):
    os.makedirs(output_path)

sigPlt.plotSBS(
    matrix_path + "breast_cancer_samples_example.SBS96.all",
    output_path,
    "BRCA_example",
    "96",
)
sigPlt.plotSBS(
    matrix_path + "breast_cancer_samples_example.SBS192.all",
    output_path,
    "BRCA_example",
    "192",
)

sigPlt.plotCNV(
    matrix_path + "breast_cancer_samples_example.CNV48.all",
    output_path,
    "BRCA_example",
    "pdf",
    percentage=False,
    aggregate=False,
)  # plotting of CNV counts
sigPlt.plotSV(
    matrix_path + "breast_cancer_samples_example.RS32.all",
    output_path,
    "BRCA_example",
    "pdf",
    percentage=False,
    aggregate=False,
)  # plotting of SV counts
