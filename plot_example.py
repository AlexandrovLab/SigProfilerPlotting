import sigProfilerPlotting as sigPlt
import os

matrix_path = "BRCA_plot/matrix/"
output_path = "BRCA_plot/plots/"

sigPlt.plot96(matrix_path + "BRCA_plot.SBS96.all", output_path, False, "BRCA_plot", False)
sigPlt.plot192(matrix_path + "BRCA_plot.SBS192.all", output_path, False, "BRCA_plot", False)
sigPlt.plotDINUC(matrix_path + "BRCA_plot.DBS78.all", output_path, False, "BRCA_plot", False)
sigPlt.plotINDEL(matrix_path + "BRCA_plot.DBS94.all", output_path, False, "BRCA_plot", False)