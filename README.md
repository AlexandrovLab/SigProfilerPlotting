# SigProfilerPlotting
SigProfilerPlotting provides a standard tool for displaying all types of mutational signatures as well as all types of mutational patterns in cancer genomes. The tool seamlessly integrates with other SigProfiler tools.

**INTRODUCTION**

The purpose of this document is to provide a guide for using the SigProfilerPlotting framework and associated functions/tools to visualize the output from SigProfilerExtraction and SigProfilerSimulator. 

**PREREQUISITES**

The framework is written in PYTHON, however, it also requires the following software with the given versions (or newer):

  * PYTHON          version 3.4 or newer
  * matplotlib [MODULE]        any version
  * SigProfilerMatrixGenerator (recommended)

**QUICK START GUIDE**

This section will guide you through the minimum steps required to create mutational matrices:
 1. Move into this repo after downloading. 
 2. Install using the following command:
 ```
 pip install .
 ```
 3.	After succesffully installing, this package can now be imported at the top of python scripts as follows:
 ```
 import sigProfilerPlotting as sigPlt 
 ```
 The available functions are listed below. 

 4. The final plots are saved into the *plots/* folder. 

**AVAILABLE FUNCTIONS**

```
import sigProfilerPlotting as sigPlt

sigPlt.plot96(matrix_path, output_path, signature, project, percentage=False)
sigPlt.plot192(matrix_path, output_path, signature, project, percentage=False)
sigPlt.plotDINUC(matrix_path, output_path, signature, project, percentage=False)
sigPlt.plotINDEL(matrix_path, output_path, signature, project, percentage=False)
```
matrix_path -> path to the mutational matrix of interest

output_path -> desired output path

signature -> Boolean: Plot based upon signature exposures

project -> name of unique sample set

percentage -> Boolean: plot the mutational matrix as percentages of the sample's total mutation count. Default is False

**EXAMPLE**

This package comes with an example test for each plot type. Run the script plot_example.py from within the downloaded repo after installation:
```
python3 plot_example.py
```
This example will create plots for each context for each of the included four samples. These plots will be saved within the BRCA_plot/plots/ folder.


**COPYRIGHT**

This software and its documentation are copyright 2018 as a part of the sigProfiler project. The sigProfilerMatrixGenerator framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Erik Bergstrom at ebergstr@eng.ucsd.edu