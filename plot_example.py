import sigProfilerPlotting as sigPlt
import os

matrix_path = "input/examples/samples/"
output_path = "output/examples/samples"

sigPlt.plot96(matrix_path + "breast_cancer_samples_example.SBS96.all", output_path, False, "BRCA_example", False)
sigPlt.plot192(matrix_path + "breast_cancer_samples_example.SBS192.all", output_path, False, "BRCA_example", False)
sigPlt.plotDINUC(matrix_path + "breast_cancer_samples_example.DBS78.all", output_path, False, "BRCA_example", False)
sigPlt.plotINDEL(matrix_path + "breast_cancer_samples_example.ID94.all", output_path, False, "BRCA_example", False)