import filecmp
import pytest
from pdf2image import convert_from_path
from PIL import Image
from PIL import ImageChops
import sigProfilerPlotting as sigPlt
import os


SPP_PATH = sigPlt.__path__[0]
SPP_SBS = os.path.join(SPP_PATH, "input/SBS/")
SPP_DBS = os.path.join(SPP_PATH, "input/DBS/")
SPP_ID = os.path.join(SPP_PATH, "input/ID/")
SPP_standard = os.path.join(SPP_PATH, "standard/")

#################
##### dumb ######
#################
# THIS SHOULD FAIL
def test_dumb_images():
    test_path = SPP_standard + "SBS_288_plots_SBS288.pdf"
    standard_path = SPP_standard + "SBS_96_plots_SBS96.pdf"

    test = convert_from_path(test_path)
    standard = convert_from_path(standard_path)

    for page_num, (test, standard) in enumerate(zip(test, standard), start=1):
        assert test.mode == standard.mode, f"Different image modes on page {page_num}."
        assert test.size == standard.size, f"Different image sizes on page {page_num}."
        assert (
            test.tobytes() == standard.tobytes()
        ), f"Images differ on page {page_num}."


def test_dumb_x_axis_images():
    test_path = SPP_standard + "SBS_288_plots_SBS288.pdf"
    standard_path = SPP_standard + "SBS_96_plots_SBS96.pdf"

    test = convert_from_path(test_path)[0]
    standard = convert_from_path(standard_path)[0]

    # Crop the images to focus on the x-axis
    x_axis_test = test.crop(
        (0, 0, test.width, 10)
    )  # Adjust the y-axis crop based on your requirement
    x_axis_standard = standard.crop(
        (0, 0, standard.width, 10)
    )  # Adjust the y-axis crop based on your requirement

    # Assert the images have the same x-axis
    assert (
        ImageChops.difference(x_axis_test, x_axis_standard).getbbox() is None
    ), f"x-axis differs"


def test_dumb_y_axis_images():
    test_path = SPP_standard + "SBS_288_plots_SBS288.pdf"
    standard_path = SPP_standard + "SBS_96_plots_SBS96.pdf"

    test = convert_from_path(test_path)[0]
    standard = convert_from_path(standard_path)[0]

    # Crop the images to focus on the y-axis
    y_axis_test = test.crop(
        (0, 0, 10, test.height)
    )  # Adjust the x-axis crop based on your requirement
    y_axis_standard = standard.crop(
        (0, 0, 10, standard.height)
    )  # Adjust the x-axis crop based on your requirement

    # Assert the images have the same y-axis
    assert (
        ImageChops.difference(y_axis_test, y_axis_standard).getbbox() is None
    ), f"y-axis differs"


def test_dumb_content_images():
    test_path = SPP_standard + "SBS_288_plots_SBS288.pdf"
    standard_path = SPP_standard + "SBS_96_plots_SBS96.pdf"

    test = convert_from_path(test_path)[0]
    standard = convert_from_path(standard_path)[0]

    # Assert the images have the same content
    assert (
        ImageChops.difference(test, standard).getbbox() is None
    ), f"Image content differs"


#################
##### SBS96 #####
#################
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
