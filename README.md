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

This section will guide you through the minimum steps required to plot mutational matrices:
 1. Before you download and install sigProfilerPlotting, first create a GitHub account at the following link: https://help.github.com/articles/signing-up-for-a-new-github-account/
 2. Configure your GitHub account to use SSH using the following link: https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/
 3. To download the current version of the tool to your local computer or server you need to clone the repository using the following command:
 ```
 git clone ssh://git@github.com/AlexandrovLab/SigProfilerPlotting 
 ```
 4. Using a terminal or command line on your computer/server, enter the sigProfilerPlotting directory
 ```
 cd sigProfilerPlotting/
 ```
 5. Install the tool using the following command:
 ```
 pip install .
 ```
 6.	After a succesfful installation, this package can now be imported at the top of python scripts as follows:
 ```
 import sigProfilerPlotting as sigPlt 
 ```
 The available functions are listed below. 

 7. The final plots are saved into the user-provided output folder. 

**ALTERNATE INSTALLATION METHOD**

SigProfilerPlotting can also be installed from pypi. Use the following command:
```
pip install sigProfilerPlotting
```

This tool can also be installed from the Anaconda Cloud as a conda package. Use the following command:
```
conda install -c ebergstr sigprofilerplotting
```

**AVAILABLE FUNCTIONS**

```
import sigProfilerPlotting as sigPlt

sigPlt.plotSBS(matrix_path, output_path, project, plot_type, percentage=False)
sigPlt.plotDBS(matrix_path, output_path, project, plot_type, percentage=False)
sigPlt.plotID(matrix_path, output_path, project, plot_type, percentage=False)
```
matrix_path -> path to the mutational matrix of interest

output_path -> desired output path

project -> name of unique sample set

plot_type -> context of the mutational matrix (96, 192, 78, 94, etc.)

percentage -> Boolean: plot the mutational matrix as percentages of the sample's total mutation count. Default is False

To create a sample portrait, ensure that you have a matrix for all required contexts (SBS-6, SBS-24, SBS-96, SBS-384, SBS-1536, SBS-6144, DBS-78, DBS-312, ID-83, ID-28, ID-96)

```
from sigProfilerPlotting import sample_portrait as sP
sP.samplePortrait(sample_matrices_path, output_path, project, percentage=False)
```
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