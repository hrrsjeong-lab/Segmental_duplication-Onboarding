"""
step00.py: Basement for everything
"""
import sys
import typing
import matplotlib.patches
import matplotlib.patheffects
import matplotlib.transforms
import numpy
import tqdm

tmpfs = "/tmpfs"
cpus_error_message = "CPUs must be positive!!"
default_error_message = "Something went wrong!!"

matplotlib_parameters = {"font.size": 50, "axes.labelsize": 50, "axes.titlesize": 60, "figure.titlesize": 60, "xtick.labelsize": 40, "ytick.labelsize": 40, "legend.fontsize": 30, "legend.title_fontsize": 30, "figure.dpi": 300, "text.color": "black", "font.family": "sans-serif", "pdf.fonttype": 42, "ps.fonttype": 42, "pdf.compression": 9}

epsilon = sys.float_info.epsilon
float_minimum = sys.float_info.min

arrowprops = {"arrowstyle": "-", "color": "silver", "linewidth": 0.5}
path_effects = [matplotlib.patheffects.withStroke(linewidth=5, foreground="white")]

sd_bed_columns = ["Chromosome 1", "Start 1", "End 1", "Name", "fakeScore", "Strand 1", "start1.1", "end1.1", "color", "Chromosome 2", "Start 2", "End 2", "Score", "Strand 2", "Maximum length", "Aligned length", "indel_a", "indel_b", "Aligned bases", "Match bases", "Mismatch bases", "Transition bases", "Transversions", "Match fraction", "Match indel fraction", "jck", "k2K", "aln_gaps", "uppercaseA", "uppercaseB", "uppercaseMatches", "aln_matches", "aln_mismatches", "aln_gap_bases", "filter_score", "sat_bases", "unique_id", "original", "telo", "peri", "acro", "telo2", "peri2", "acro2"]


def check_suffix(filename: str, suffixes: typing.Set[str]) -> None:
    filename = filename.lower()
    for suffix in suffixes:
        if filename.endswith(suffix):
            return
    raise ValueError(f"{filename} must be ended with one of {sorted(suffixes)}!!")


def check_suffixes(filenames: typing.List[str], suffixes: typing.Set[str]) -> None:
    for filename in tqdm.tqdm(filenames):
        check_suffix(filename, suffixes)


def confidence_ellipse(x: typing.List[float], y: typing.List[float], ax, n_std: float = 2.0, facecolor: str = "none", **kwargs) -> matplotlib.patches.Patch:
    if len(x) != len(y):
        raise ValueError("x and y must be the same size!!")

    if len(x) < 3:
        return matplotlib.patches.Ellipse((0, 0), width=0, height=0, facecolor=facecolor, **kwargs)

    cov = numpy.cov(x, y)
    pearson = cov[0, 1] / numpy.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = numpy.sqrt(1 + pearson)
    ell_radius_y = numpy.sqrt(1 - pearson)
    ellipse = matplotlib.patches.Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, facecolor=facecolor, **kwargs)

    scale_x = numpy.sqrt(cov[0, 0]) * n_std
    mean_x = numpy.mean(x)

    scale_y = numpy.sqrt(cov[1, 1]) * n_std
    mean_y = numpy.mean(y)

    transf = matplotlib.transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def pvalue_format(p_value: float) -> str:
    thresholds = [1e-4, 1e-3, 1e-2, 5e-2]

    if not (0.0 <= p_value <= 1.0):
        raise ValueError(f"p={p_value} is not a valid p-value!!")

    if p_value > thresholds[-1]:
        return f"p={p_value:.3f}"
    elif p_value < thresholds[0]:
        return f"p={p_value:.1e}"

    for threshold in thresholds:
        if p_value < threshold:
            return f"p<{p_value:.3e}"

    raise ValueError(default_error_message)


def format_filename(s: str) -> str:
    return s.title().replace(" ", "")
