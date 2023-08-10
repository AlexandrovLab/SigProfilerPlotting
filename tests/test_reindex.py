import os

from PIL import Image, ImageChops

import sigProfilerPlotting as sigPlt

current_script_path = os.path.abspath(__file__)
SPP_PATH = os.path.dirname(current_script_path)

SPP_SBS = os.path.join(SPP_PATH, "input/SBS/")
SPP_DBS = os.path.join(SPP_PATH, "input/DBS/")
SPP_ID = os.path.join(SPP_PATH, "input/ID/")
SPP_STANDARD_PNG = os.path.join(SPP_PATH, "standard_png/")


def image_difference(img1, img2):
    """
    Calculate the difference between two images using ImageChops.

    Parameters:
    - img1_path: Path to the first image.
    - img2_path: Path to the second image.

    Returns:
    - relative_difference: The relative difference between the two images.
    """

    img1_opened_in_func = False
    img2_opened_in_func = False

    # Open the images if they are not already opened
    if isinstance(img1, str):
        img1 = Image.open(img1)
        img1_opened_in_func = True
    if isinstance(img2, str):
        img2 = Image.open(img2)
        img2_opened_in_func = True

    # Convert images to grayscale
    img1 = img1.convert("L")
    img2 = img2.convert("L")

    # Difference between the two images
    diff = ImageChops.difference(img1, img2)

    # Calculate the absolute difference
    total_difference = sum(abs(p) for p in diff.getdata())
    width, height = img1.size
    num_channels = 1
    max_difference = width * height * num_channels * 255

    # Calculate the relative difference
    relative_difference = total_difference / max_difference

    # Close the images if they were opened
    if img1_opened_in_func:
        img1.close()
    if img2_opened_in_func:
        img2.close()

    return relative_difference


#################
##### SBS96 #####
#################
# test full image using unordered input
def test_SBS96_unordered_images():
    # creat the output directory path
    output_directory = os.path.join(SPP_SBS, "output", "96_full_image", "")
    os.makedirs(output_directory, exist_ok=True)

    sigPlt.plotSBS(
        SPP_SBS + "unordered/example.SBS96.all",
        output_directory,
        "test_unordered",
        "96",
        savefig_format="png",
    )

    test_path = os.path.join(output_directory, "SBS_96_plots_Random.png")
    standard_path = os.path.join(SPP_STANDARD_PNG, "SBS_96_plots_Random.png")

    # compare the two images
    relative_difference = image_difference(test_path, standard_path)

    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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
    x_axis_test = test.crop((170, 905, 6000, 2000))

    relative_difference = image_difference(x_axis_test, standard_path)
    test.close()

    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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

    # Crop the images to focus on the bars
    test = Image.open(test_path)
    bars_test = test.crop((170, 300, 4340, 910))

    # Compare the images
    relative_difference = image_difference(bars_test, standard_path)
    test.close()

    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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

    relative_difference = image_difference(test_path, standard_path)

    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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

    # Crop the images to focus on the x-axis
    test = Image.open(test_path)
    x_axis_test = test.crop((170, 905, 3250, 2000))

    relative_difference = image_difference(x_axis_test, standard_path)
    test.close()

    # Compare the images
    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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

    # Crop the images to focus on the bars
    test = Image.open(test_path)
    bars_test = test.crop((170, 300, 3250, 910))

    relative_difference = image_difference(bars_test, standard_path)
    test.close()

    # Compare the images
    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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

    relative_difference = image_difference(test_path, standard_path)

    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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

    # Crop the images to focus on the x-axis
    test = Image.open(test_path)
    x_axis_test = test.crop((170, 905, 4350, 990))

    relative_difference = image_difference(x_axis_test, standard_path)
    test.close()

    # Compare the images
    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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

    # Crop the images to focus on the bars
    test = Image.open(test_path)
    bars_test = test.crop((170, 280, 4350, 910))

    # Compare the images
    relative_difference = image_difference(bars_test, standard_path)
    test.close()

    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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

    relative_difference = image_difference(test_path, standard_path)

    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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

    # Crop the images to focus on the x-axis
    test = Image.open(test_path)
    x_axis_test = test.crop((170, 1063, 4250, 1100))

    # Compare the images
    relative_difference = image_difference(x_axis_test, standard_path)
    test.close()

    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"


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

    # Crop the images to focus on the bars
    test = Image.open(test_path)
    bars_test = test.crop((190, 350, 4250, 1003))

    relative_difference = image_difference(bars_test, standard_path)
    test.close()

    assert (
        relative_difference == 0
    ), f"relative difference between unordered and standard is {relative_difference}"
