import sigProfilerPlotting as sigPlt
import os

matrix_path = "input/examples/samples/"
output_path = "output/examples/samples"

cnv_matrix_path = ""
sv_Matrix_path = ""

#matrix_path = ''
#output_path = ''



sigPlt.plotSBS(matrix_path + "breast_cancer_samples_example.SBS96.all", output_path,"BRCA_example", "96")
sigPlt.plotSBS(matrix_path + "breast_cancer_samples_example.SBS192.all", output_path, "BRCA_example", "192")

sigPlt.plotCNV(matrix_path + "breast_cancer_samples_example.CNV48.all", output_path, "BRCA_example", "pdf", percentage=False, aggregate=False)
sigPlt.plotSV(matrix_path + "breast_cancer_samples_example.RS32.all", output_path, "BRCA_example", "pdf", percentage=False, aggregate=False)

#sigPlt.plotDBS(matrix_path + "breast_cancer_samples_example.DBS78.all", output_path, "BRCA_example", "78")
#sigPlt.plotID(matrix_path + "breast_cancer_samples_example.ID94.all", output_path,  "BRCA_example", "94")
