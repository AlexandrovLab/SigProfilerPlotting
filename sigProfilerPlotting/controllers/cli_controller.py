import argparse
from typing import List
import sigProfilerPlotting as sigPlt


# Common parser setup for shared arguments
def common_plotting_arguments(parser):
    parser.add_argument("matrix_path", help="The path to the input matrix file.")
    parser.add_argument(
        "output_path", help="The directory where the plots will be saved."
    )
    parser.add_argument("project", help="The name of the project.")
    parser.add_argument("plot_type", help="The type of plot to generate.")
    parser.add_argument(
        "--percentage", action="store_true", help="Display percentages in the plot."
    )
    parser.add_argument(
        "--custom_text_upper", help="Custom text to display at the top of the plot."
    )
    parser.add_argument(
        "--custom_text_middle", help="Custom text to display in the middle of the plot."
    )
    parser.add_argument(
        "--custom_text_bottom", help="Custom text to display at the bottom of the plot."
    )
    parser.add_argument(
        "--savefig_format",
        default="pdf",
        choices=["pdf", "png", "pil_image"],
        help="The file format for saving the plot.",
    )
    parser.add_argument("--volume", help="Specify a volume for Docker container usage.")
    parser.add_argument(
        "--dpi",
        type=int,
        default=100,
        help="The resolution of the plot in dots per inch.",
    )


def parse_arguments_sbs(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="SigProfilerPlotting plotDBS", description="Generate DBS plots."
    )
    common_plotting_arguments(parser)
    return parser.parse_args(args)


def parse_arguments_id(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="SigProfilerPlotting plotID", description="Generate ID plots."
    )
    common_plotting_arguments(parser)
    return parser.parse_args(args)


def parse_arguments_dbs(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="SigProfilerPlotting plotDBS", description="Generate DBS plots."
    )
    common_plotting_arguments(parser)
    return parser.parse_args(args)


def parse_arguments_sv(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="SigProfilerPlotting plotSV", description="Generate SV plots."
    )
    parser.add_argument("matrix_path", help="The path to the input matrix file.")
    parser.add_argument(
        "output_path", help="The directory where the plots will be saved."
    )
    parser.add_argument("project", help="The name of the project.")
    parser.add_argument(
        "--percentage", action="store_true", help="Display percentages in the plot."
    )
    parser.add_argument(
        "--aggregate", action="store_true", help="Aggregate all samples."
    )
    parser.add_argument(
        "--savefig_format",
        default="pdf",
        choices=["pdf", "png", "pil_image"],
        help="The file format for saving the plot.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=100,
        help="The resolution of the plot in dots per inch.",
    )
    return parser.parse_args(args)


def parse_arguments_cnv(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="SigProfilerPlotting plotCNV", description="Generate CNV plots."
    )
    parser.add_argument("matrix_path", help="The path to the input matrix file.")
    parser.add_argument(
        "output_path", help="The directory where the plots will be saved."
    )
    parser.add_argument("project", help="The name of the project.")
    parser.add_argument(
        "--percentage", action="store_true", help="Display percentages in the plot."
    )
    parser.add_argument(
        "--aggregate", action="store_true", help="Aggregate data for the plot."
    )
    parser.add_argument(
        "--read_from_file",
        action="store_true",
        default=True,
        help="Read data from a file for the plot.",
    )
    parser.add_argument(
        "--savefig_format",
        default="pdf",
        choices=["pdf", "png", "pil_image"],
        help="The file format for saving the plot.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=100,
        help="The resolution of the plot in dots per inch.",
    )
    return parser.parse_args(args)


def dispatch_plot_sbs(parsed_args: argparse.Namespace) -> None:
    sigPlt.plotSBS(
        matrix_path=parsed_args.matrix_path,
        output_path=parsed_args.output_path,
        project=parsed_args.project,
        plot_type=parsed_args.plot_type,
        percentage=parsed_args.percentage,
        custom_text_upper=parsed_args.custom_text_upper,
        custom_text_middle=parsed_args.custom_text_middle,
        custom_text_bottom=parsed_args.custom_text_bottom,
        savefig_format=parsed_args.savefig_format,
        volume=parsed_args.volume,
        dpi=parsed_args.dpi,
    )


def dispatch_plot_id(parsed_args: argparse.Namespace) -> None:
    sigPlt.plotID(
        matrix_path=parsed_args.matrix_path,
        output_path=parsed_args.output_path,
        project=parsed_args.project,
        plot_type=parsed_args.plot_type,
        percentage=parsed_args.percentage,
        custom_text_upper=parsed_args.custom_text_upper,
        custom_text_middle=parsed_args.custom_text_middle,
        custom_text_bottom=parsed_args.custom_text_bottom,
        savefig_format=parsed_args.savefig_format,
        volume=parsed_args.volume,
        dpi=parsed_args.dpi,
    )


def dispatch_plot_dbs(parsed_args: argparse.Namespace) -> None:
    sigPlt.plotDBS(
        matrix_path=parsed_args.matrix_path,
        output_path=parsed_args.output_path,
        project=parsed_args.project,
        plot_type=parsed_args.plot_type,
        percentage=parsed_args.percentage,
        custom_text_upper=parsed_args.custom_text_upper,
        custom_text_middle=parsed_args.custom_text_middle,
        custom_text_bottom=parsed_args.custom_text_bottom,
        savefig_format=parsed_args.savefig_format,
        volume=parsed_args.volume,
        dpi=parsed_args.dpi,
    )


def dispatch_plot_sv(parsed_args: argparse.Namespace) -> None:
    sigPlt.plotSV(
        matrix_path=parsed_args.matrix_path,
        output_path=parsed_args.output_path,
        project=parsed_args.project,
        percentage=parsed_args.percentage,
        aggregate=parsed_args.aggregate,
        savefig_format=parsed_args.savefig_format,
        dpi=parsed_args.dpi,
    )


def dispatch_plot_cnv(parsed_args: argparse.Namespace) -> None:
    sigPlt.plotCNV(
        matrix_path=parsed_args.matrix_path,
        output_path=parsed_args.output_path,
        project=parsed_args.project,
        percentage=parsed_args.percentage,
        aggregate=parsed_args.aggregate,
        read_from_file=parsed_args.read_from_file,
        savefig_format=parsed_args.savefig_format,
        dpi=parsed_args.dpi,
    )


class CliController:
    def dispatch(self, user_args: List[str]):
        if "plotSBS" in user_args:
            parsed_args = parse_arguments_sbs(user_args[1:])
            dispatch_plot_sbs(parsed_args)
        elif "plotID" in user_args:
            parsed_args = parse_arguments_id(user_args[1:])
            dispatch_plot_id(parsed_args)
        elif "plotDBS" in user_args:
            parsed_args = parse_arguments_dbs(user_args[1:])
            dispatch_plot_dbs(parsed_args)
        elif "plotSV" in user_args:
            parsed_args = parse_arguments_sv(user_args[1:])
            dispatch_plot_sv(parsed_args)
        elif "plotCNV" in user_args:
            parsed_args = parse_arguments_cnv(user_args[1:])
            dispatch_plot_cnv(parsed_args)
        else:
            print(
                "Unknown command. Available commands: plotSBS, plotID, plotDBS, plotSV, plotCNV."
            )


if __name__ == "__main__":
    import sys

    controller = CliController()
    controller.dispatch(sys.argv)
