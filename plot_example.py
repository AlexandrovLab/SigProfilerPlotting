import sigProfilerPlotting as sigPlt
import os

matrix_path = "input/examples/"
output_path = "output/examples/"

#matrix_path = ''
#output_path = ''

sigPlt.plotSBS(matrix_path + "breast_cancer_samples_example.SBS96.all", output_path,"BRCA_example", "96")
sigPlt.plotSBS(matrix_path + "breast_cancer_samples_example.SBS192.all", output_path, "BRCA_example", "192")
#sigPlt.plotDBS(matrix_path + "breast_cancer_samples_example.DBS78.all", output_path, "BRCA_example", "78")
#sigPlt.plotID(matrix_path + "breast_cancer_samples_example.ID94.all", output_path,  "BRCA_example", "94")
