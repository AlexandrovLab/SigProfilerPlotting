import filecmp
import pytest
from PIL import Image
from PIL import ImageChops
import sigProfilerPlotting as sigPlt
import os

# To run locally:`/Users/tingyang/opt/anaconda3/envs/u56/bin/python -m pytest tests`

current_script_path = os.path.abspath(__file__)
SPP_PATH = os.path.dirname(current_script_path)

SPP_SBS = os.path.join(SPP_PATH, "input/SBS/")
SPP_DBS = os.path.join(SPP_PATH, "input/DBS/")
SPP_ID = os.path.join(SPP_PATH, "input/ID/")
SPP_STANDARD_PNG = os.path.join(SPP_PATH, "standard_png/")

#################
##### SBS96 #####
#################
# test full image using unordered input
def test_SBS96_unordered_images():
    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS96.all",
        SPP_SBS + "output/96_full_image/",
        "test_unordered",
        "96",
        savefig_format="png",
    )

    test_path = SPP_SBS + "output/96_full_image/SBS_96_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "SBS_96_plots_Random.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    assert test.tobytes() == standard.tobytes(), f"image differs"


# test x-axis using unordered input
def test_SBS96_x_axis_images():
    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS96.all",
        SPP_SBS + "output/96_x_axis/",
        "test_unordered",
        "96",
        savefig_format="png",
    )

    test_path = SPP_SBS + "output/96_x_axis/SBS_96_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "SBS_96_plots_xaxis.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the x-axis
    x_axis_test = test.crop((170, 905, 6000, 2000))

    # Compare the images
    assert x_axis_test.tobytes() == standard.tobytes(), f"x-axis differs"


# test bars using unordered input
def test_SBS96_bars_images():
    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS96.all",
        SPP_SBS + "output/96_bars/",
        "test_unordered",
        "96",
        savefig_format="png",
    )

    test_path = SPP_SBS + "output/96_bars/SBS_96_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "SBS_96_plots_bars.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the bars
    bars_test = test.crop((170, 300, 4340, 910))

    # Compare the images
    assert bars_test.tobytes() == standard.tobytes(), f"bars differs"


##################
##### SBS288 #####
##################
# test full image using unordered input
def test_SBS288_unordered_images():
    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS288.all",
        SPP_SBS + "output/288_full_image/",
        "test_unordered",
        "288",
        savefig_format="png",
    )

    test_path = SPP_SBS + "output/288_full_image/SBS_288_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "SBS_288_plots_Random.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    assert test.tobytes() == standard.tobytes(), f"image differs"


# test x-axis using unordered input
def test_SBS288_x_axis_images():
    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS288.all",
        SPP_SBS + "output/288_x_axis/",
        "test_unordered",
        "288",
        savefig_format="png",
    )

    test_path = SPP_SBS + "output/288_x_axis/SBS_288_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "SBS_288_plots_xaxis.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the x-axis
    x_axis_test = test.crop((170, 905, 3250, 2000))

    # Compare the images
    assert x_axis_test.tobytes() == standard.tobytes(), f"x-axis differs"


# test bars using unordered input
def test_SBS288_bars_images():
    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS288.all",
        SPP_SBS + "output/288_bars/",
        "test_unordered",
        "288",
        savefig_format="png",
    )

    test_path = SPP_SBS + "output/288_bars/SBS_288_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "SBS_288_plots_bars.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the bars
    bars_test = test.crop((170, 300, 3250, 910))

    # Compare the images
    assert bars_test.tobytes() == standard.tobytes(), f"bars differs"


#################
##### DBS78 #####
#################
# test full image using unordered input
def test_DBS78_unordered_images():
    sigPlt.plotDBS(
        SPP_DBS + "unordered/example.DBS78.all",
        SPP_DBS + "output/78_full_image/",
        "test_unordered",
        "78",
        savefig_format="png",
    )

    test_path = SPP_DBS + "output/78_full_image/DBS_78_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "DBS_78_plots_Random.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    assert test.tobytes() == standard.tobytes(), f"image differs"


# test x-axis using unordered input
def test_DBS78_x_axis_images():
    sigPlt.plotDBS(
        SPP_DBS + "unordered/example.DBS78.all",
        SPP_DBS + "output/78_x_axis/",
        "test_unordered",
        "78",
        savefig_format="png",
    )

    test_path = SPP_DBS + "output/78_x_axis/DBS_78_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "DBS_78_plots_xaxis.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the x-axis
    x_axis_test = test.crop((170, 905, 4350, 990))

    # Compare the images
    assert x_axis_test.tobytes() == standard.tobytes(), f"x-axis differs"


# test bars using unordered input
def test_DBS78_bars_images():
    sigPlt.plotDBS(
        SPP_DBS + "unordered/example.DBS78.all",
        SPP_DBS + "output/78_bars/",
        "test_unordered",
        "78",
        savefig_format="png",
    )

    test_path = SPP_DBS + "output/78_bars/DBS_78_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "DBS_78_plots_bars.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the bars
    bars_test = test.crop((170, 280, 4350, 910))

    # Compare the images
    assert bars_test.tobytes() == standard.tobytes(), f"bars differs"


################
##### ID83 #####
################
# test full image using unordered input
def test_ID83_unordered_images():
    sigPlt.plotID(
        SPP_ID + "unordered/example.ID83.all",
        SPP_ID + "output/83_full_image/",
        "test_unordered",
        "83",
        savefig_format="png",
    )

    test_path = SPP_ID + "output/83_full_image/ID_83_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "ID_83_plots_Random.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    assert test.tobytes() == standard.tobytes(), f"image differs"


# test x-axis using unordered input
def test_ID83_x_axis_images():
    sigPlt.plotID(
        SPP_ID + "unordered/example.ID83.all",
        SPP_ID + "output/83_x_axis/",
        "test_unordered",
        "83",
        savefig_format="png",
    )

    test_path = SPP_ID + "output/83_x_axis/ID_83_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "ID_83_plots_xaxis.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the x-axis
    x_axis_test = test.crop((170, 1063, 4250, 1100))

    # Compare the images
    assert x_axis_test.tobytes() == standard.tobytes(), f"x-axis differs"


# test bars using unordered input
def test_ID83_bars_images():
    sigPlt.plotID(
        SPP_ID + "unordered/example.ID83.all",
        SPP_ID + "output/83_bars/",
        "test_unordered",
        "83",
        savefig_format="png",
    )

    test_path = SPP_ID + "output/83_bars/ID_83_plots_Random.png"
    standard_path = SPP_STANDARD_PNG + "ID_83_plots_bars.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the bars
    bars_test = test.crop((190, 350, 4250, 1003))

    # Compare the images
    assert bars_test.tobytes() == standard.tobytes(), f"bars differs"
