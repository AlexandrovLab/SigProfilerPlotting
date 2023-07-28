import filecmp
import pytest
from pdf2image import convert_from_path
from PIL import Image
from PIL import ImageChops
import sigProfilerPlotting as sigPlt
import os

current_script_path = os.path.abspath(__file__)
SPP_PATH = os.path.dirname(current_script_path)

# SPP_PATH = sigPlt.__path__[0]
SPP_SBS = os.path.join(SPP_PATH, "input/SBS/")
SPP_DBS = os.path.join(SPP_PATH, "input/DBS/")
SPP_ID = os.path.join(SPP_PATH, "input/ID/")
SPP_standard = os.path.join(SPP_PATH, "standard/")
SPP_standard_png = os.path.join(SPP_PATH, "standard_png/")

#################
##### SBS96 #####
#################
# test full image using ordered input
def test_SBS96_images():
    sigPlt.plotSBS(
        SPP_SBS + "ordered/example.SBS96.all", SPP_SBS + "output/", "test_ordered", "96"
    )

    test_path = SPP_SBS + "output/SBS_96_plots_test_ordered.pdf"
    standard_path = SPP_standard + "SBS_96_plots_SBS96.pdf"

    test = convert_from_path(test_path)
    standard = convert_from_path(standard_path)

    for page_num, (test, standard) in enumerate(zip(test, standard), start=1):
        assert test.mode == standard.mode, f"Different image modes on page {page_num}."
        assert test.size == standard.size, f"Different image sizes on page {page_num}."
        assert (
            test.tobytes() == standard.tobytes()
        ), f"Images differ on page {page_num}."


# test full image using unordered input
def test_SBS96_unordered_images():
    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS96.all",
        SPP_SBS + "output/",
        "test_unordered",
        "96",
    )

    test_path = SPP_SBS + "output/SBS_96_plots_test_unordered.pdf"
    standard_path = SPP_standard + "SBS_96_plots_SBS96.pdf"

    test = convert_from_path(test_path)
    standard = convert_from_path(standard_path)

    for page_num, (test, standard) in enumerate(zip(test, standard), start=1):
        assert test.mode == standard.mode, f"Different image modes on page {page_num}."
        assert test.size == standard.size, f"Different image sizes on page {page_num}."
        assert (
            test.tobytes() == standard.tobytes()
        ), f"Images differ on page {page_num}."


# test x-axis using unordered input
def test_SBS96_x_axis_images():
    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS96.all",
        SPP_SBS + "output/",
        "test_unordered",
        "96",
        savefig_format="png",
    )

    test_path = SPP_SBS + "output/SBS_96_plots_Random.png"
    standard_path = SPP_standard_png + "SBS_96_plots_xaxis.png"

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
        SPP_SBS + "output/",
        "test_unordered",
        "96",
        savefig_format="png",
    )

    test_path = SPP_SBS + "output/SBS_96_plots_Random.png"
    standard_path = SPP_standard_png + "SBS_96_plots_bars.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the bars
    bars_test = test.crop((170, 300, 4340, 910))

    # Compare the images
    assert bars_test.tobytes() == standard.tobytes(), f"bars differs"


##################
##### SBS288 #####
##################
def test_SBS288_images():
    sigPlt.plotSBS(
        SPP_SBS + "ordered/example.SBS288.all",
        SPP_SBS + "output/",
        "test_ordered",
        "288",
    )

    test_path = SPP_SBS + "output/SBS_288_plots_test_ordered.pdf"
    standard_path = SPP_standard + "SBS_288_plots_SBS288.pdf"

    test = convert_from_path(test_path)
    standard = convert_from_path(standard_path)

    for page_num, (test, standard) in enumerate(zip(test, standard), start=1):
        assert test.mode == standard.mode, f"Different image modes on page {page_num}."
        assert test.size == standard.size, f"Different image sizes on page {page_num}."
        assert (
            test.tobytes() == standard.tobytes()
        ), f"Images differ on page {page_num}."


def test_SBS288_unordered_images():
    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS288.all",
        SPP_SBS + "output/",
        "test_unordered",
        "288",
    )

    test_path = SPP_SBS + "output/SBS_288_plots_test_unordered.pdf"
    standard_path = SPP_standard + "SBS_288_plots_SBS288.pdf"

    test = convert_from_path(test_path)
    standard = convert_from_path(standard_path)

    for page_num, (test, standard) in enumerate(zip(test, standard), start=1):
        assert test.mode == standard.mode, f"Different image modes on page {page_num}."
        assert test.size == standard.size, f"Different image sizes on page {page_num}."
        assert (
            test.tobytes() == standard.tobytes()
        ), f"Images differ on page {page_num}."


# test x-axis using unordered input
def test_SBS288_x_axis_images():
    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS288.all",
        SPP_SBS + "output/",
        "test_unordered",
        "288",
        savefig_format="png",
    )

    test_path = SPP_SBS + "output/SBS_288_plots_Random.png"
    standard_path = SPP_standard_png + "SBS_288_plots_xaxis.png"

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
        SPP_SBS + "output/",
        "test_unordered",
        "288",
        savefig_format="png",
    )

    test_path = SPP_SBS + "output/SBS_288_plots_Random.png"
    standard_path = SPP_standard_png + "SBS_288_plots_bars.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the bars
    bars_test = test.crop((170, 300, 3250, 910))

    # Compare the images
    assert bars_test.tobytes() == standard.tobytes(), f"bars differs"


#################
##### DBS78 #####
#################
def test_DBS78_images():
    sigPlt.plotDBS(
        SPP_DBS + "ordered/example.DBS78.all",
        SPP_DBS + "output/",
        "test_ordered",
        "78",
    )

    test_path = SPP_DBS + "output/DBS_78_plots_test_ordered.pdf"
    standard_path = SPP_standard + "DBS_78_plots_DBS78.pdf"

    test = convert_from_path(test_path)
    standard = convert_from_path(standard_path)

    for page_num, (test, standard) in enumerate(zip(test, standard), start=1):
        assert test.mode == standard.mode, f"Different image modes on page {page_num}."
        assert test.size == standard.size, f"Different image sizes on page {page_num}."
        assert (
            test.tobytes() == standard.tobytes()
        ), f"Images differ on page {page_num}."


def test_DBS78_unordered_images():
    sigPlt.plotDBS(
        SPP_DBS + "unordered/example.DBS78.all",
        SPP_DBS + "output/",
        "test_unordered",
        "78",
    )

    test_path = SPP_DBS + "output/DBS_78_plots_test_unordered.pdf"
    standard_path = SPP_standard + "DBS_78_plots_DBS78.pdf"

    test = convert_from_path(test_path)
    standard = convert_from_path(standard_path)

    for page_num, (test, standard) in enumerate(zip(test, standard), start=1):
        assert test.mode == standard.mode, f"Different image modes on page {page_num}."
        assert test.size == standard.size, f"Different image sizes on page {page_num}."
        assert (
            test.tobytes() == standard.tobytes()
        ), f"Images differ on page {page_num}."


# test x-axis using unordered input
def test_DBS78_x_axis_images():
    sigPlt.plotDBS(
        SPP_DBS + "unordered/example.DBS78.all",
        SPP_DBS + "output/",
        "test_unordered",
        "78",
        savefig_format="png",
    )

    test_path = SPP_DBS + "output/DBS_78_plots_Random.png"
    standard_path = SPP_standard_png + "DBS_78_plots_xaxis.png"

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
        SPP_DBS + "output/",
        "test_unordered",
        "78",
        savefig_format="png",
    )

    test_path = SPP_DBS + "output/DBS_78_plots_Random.png"
    standard_path = SPP_standard_png + "DBS_78_plots_bars.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the bars
    bars_test = test.crop((170, 280, 4350, 910))

    # Compare the images
    assert bars_test.tobytes() == standard.tobytes(), f"bars differs"


################
##### ID83 #####
################
def test_ID83_images():
    sigPlt.plotID(
        SPP_ID + "ordered/example.ID83.all",
        SPP_ID + "output/",
        "test_ordered",
        "83",
    )

    test_path = SPP_ID + "output/ID_83_plots_test_ordered.pdf"
    standard_path = SPP_standard + "ID_83_plots_ID83.pdf"

    test = convert_from_path(test_path)
    standard = convert_from_path(standard_path)

    for page_num, (test, standard) in enumerate(zip(test, standard), start=1):
        assert test.mode == standard.mode, f"Different image modes on page {page_num}."
        assert test.size == standard.size, f"Different image sizes on page {page_num}."
        assert (
            test.tobytes() == standard.tobytes()
        ), f"Images differ on page {page_num}."


def test_ID83_unordered_images():
    sigPlt.plotID(
        SPP_ID + "unordered/example.ID83.all",
        SPP_ID + "output/",
        "test_unordered",
        "83",
    )

    test_path = SPP_ID + "output/ID_83_plots_test_unordered.pdf"
    standard_path = SPP_standard + "ID_83_plots_ID83.pdf"

    test = convert_from_path(test_path)
    standard = convert_from_path(standard_path)

    for page_num, (test, standard) in enumerate(zip(test, standard), start=1):
        assert test.mode == standard.mode, f"Different image modes on page {page_num}."
        assert test.size == standard.size, f"Different image sizes on page {page_num}."
        assert (
            test.tobytes() == standard.tobytes()
        ), f"Images differ on page {page_num}."


# test x-axis using unordered input
def test_ID83_x_axis_images():
    sigPlt.plotID(
        SPP_ID + "unordered/example.ID83.all",
        SPP_ID + "output/",
        "test_unordered",
        "83",
        savefig_format="png",
    )

    test_path = SPP_ID + "output/ID_83_plots_Random.png"
    standard_path = SPP_standard_png + "ID_83_plots_xaxis.png"

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
        SPP_ID + "output/",
        "test_unordered",
        "83",
        savefig_format="png",
    )

    test_path = SPP_ID + "output/ID_83_plots_Random.png"
    standard_path = SPP_standard_png + "ID_83_plots_bars.png"

    test = Image.open(test_path)
    standard = Image.open(standard_path)

    # Crop the images to focus on the bars
    bars_test = test.crop((190, 350, 4250, 1003))

    # Compare the images
    assert bars_test.tobytes() == standard.tobytes(), f"bars differs"
