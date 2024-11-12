import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import re
import datetime
import shutil
import sys
import pathlib

def set_matplotlib_style():

    textsize = 16
    markersize = 4.0
    labelsize =  12.5

    plt.rcParams.update({
        # "text.usetex": True,
        # "font.family": "serif",
        # "font.serif": ["Palatino",],
        "font.size": textsize,
        'text.latex.preamble': [r"""\usepackage{bm}"""]
    })

    mpl.rcParams.update({
        # "legend.fontsize" : "medium",
        "legend.frameon" : False,
        "legend.handletextpad" : 0.1,
        "legend.columnspacing" : 0.8,
        # "axes.labelsize" : "medium",
        # "axes.titlesize" : "medium",
        # "xtick.labelsize" : "medium",
        # "ytick.labelsize" : "medium",
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
        'xtick.minor.visible': True,
        'ytick.minor.visible': True,
    })
    colors = ["#e74c3c","#f1c40f","#16a085","#27ae60","#2980b9","#8e44ad"]

    return colors, textsize, labelsize, markersize


def script_command_to_str(argv, parser_args):
    call_args = np.array(argv[1:], dtype=object)
    match_expr = "|".join(["^-+([a-z]+[1-9]*-*)+"]+([] if not parser_args else [f"^-*{x.replace('_', '.')}" for x in vars(parser_args).keys()]))
    if call_args.size != 0:
        flags = np.vectorize(lambda x: bool(re.match(match_expr, x)))(call_args)
        special_chars = np.vectorize(lambda x: not x.isalnum())(call_args)
        select = ~flags & special_chars
        if np.count_nonzero(select):
            call_args[select] = np.vectorize(lambda x: f"'{x}'")(call_args[select])
    return " ".join([argv[0], *call_args])


def write_index_and_log(
    outpath,
    logname,
    template_dir=f"{pathlib.Path(__file__).parent}/../res/",
    yield_tables=None,
    analysis_meta_info=None,
    args={},
    nround=2,
):
    indexname = "index.php"
    shutil.copyfile(f"{template_dir}/{indexname}", f"{outpath}/index.php")
    logname = f"{outpath}/{logname}.log"

    with open(logname, "w") as logf:
        meta_info = (
            "-" * 80
            + "\n"
            + f"Script called at {datetime.datetime.now()}\n"
            + f"The command was: {script_command_to_str(sys.argv, args)}\n"
            + "-" * 80
            + "\n"
        )
        logf.write(meta_info)

        if yield_tables:
            for k, v in yield_tables.items():
                logf.write(f"Yield information for {k}\n")
                logf.write("-" * 80 + "\n")
                logf.write(str(v.round(nround)) + "\n\n")

            if (
                "Unstacked processes" in yield_tables
                and "Stacked processes" in yield_tables
            ):
                if "Data" in yield_tables["Unstacked processes"]["Process"].values:
                    unstacked = yield_tables["Unstacked processes"]
                    data_yield = unstacked[unstacked["Process"] == "Data"][
                        "Yield"
                    ].iloc[0]
                    ratio = (
                        float(
                            yield_tables["Stacked processes"]["Yield"].sum()
                            / data_yield
                        )
                        * 100
                    )
                    logf.write(f"===> Sum unstacked to data is {ratio:.2f}%")

        if analysis_meta_info:
            for k, analysis_info in analysis_meta_info.items():
                logf.write("\n" + "-" * 80 + "\n")
                logf.write(f"Meta info from input file {k}\n")
                logf.write("\n" + "-" * 80 + "\n")
                logf.write(json.dumps(analysis_info, indent=5).replace("\\n", "\n"))
        print(f"Writing file {logname}")