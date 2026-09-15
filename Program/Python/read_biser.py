import argparse
import itertools
import os
import tempfile
import typing
import zipfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input BEDPE file", type=str)
    parser.add_argument("reference", help="Reference BED(.gz) file", type=str)
    parser.add_argument("length", help="Length TXT file", type=str)
    parser.add_argument("output", help="Output ZIP file", type=str)

    args = parser.parse_args()

    step00.check_suffix(args.input, {".bedpe"})
    step00.check_suffix(args.reference, {".bed", ".bed.gz"})
    step00.check_suffix(args.length, {".txt"})
    step00.check_suffix(args.output, {".zip"})

    column_data_type_dict = {"Chromosome 1": str, "Start 1": int, "End 1": int, "Chromosome 2": str, "Start 2": int, "End 2": int, "Name": str, "Score": float, "Strand 1": str, "Strand 2": str, "Maximum length": int, "Aligned length": int, "cigar": str, "comment": str}
    type_palette = {"BISER": "tab:blue", "Reference": "tab:orange"}
    type_marker = {"BISER": "P", "Reference": "X"}

    input_data = pandas.read_csv(args.input, sep="\t", header=None).iloc[:, :len(column_data_type_dict)]
    input_data.columns = list(column_data_type_dict.keys())
    input_data = input_data.astype(column_data_type_dict).sort_values(["Chromosome 1", "Start 1"])
    print(input_data)

    input_data["Type"] = "BISER"
    print(input_data)

    chromosome_pair_set = set(input_data[["Chromosome 1", "Chromosome 2"]].itertuples(index=False, name=None))
    print("Chromosome pairs:", len(chromosome_pair_set))

    reference_data = pandas.read_csv(args.reference, sep="\t")
    reference_data.columns = step00.sd_bed_columns
    reference_data = reference_data.sort_values(["Chromosome 1", "Start 1"])
    reference_data["Type"] = "Reference"
    print(reference_data)

    chromosome_pair_set &= set(reference_data[["Chromosome 1", "Chromosome 2"]].itertuples(index=False, name=None))
    print("Chromosome pairs:", len(chromosome_pair_set))

    length_data = pandas.read_csv(args.length, sep="\t", header=None, names=["Chromosome", "Length"], index_col=0)
    chromosome_list = list(length_data.index)
    print(length_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with tempfile.TemporaryDirectory() as directory:
        figure_list: typing.List[str] = list()

        drawing_data = pandas.concat([input_data, reference_data], join="inner", ignore_index=True)

        # ECDF plots
        for x in tqdm.tqdm(["Score"], position=1):
            fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

            seaborn.ecdfplot(data=drawing_data, x=x, stat="percent", hue="Type", palette=type_palette, ax=ax)

            matplotlib.pyplot.tight_layout()

            figure_list.append(f"{directory}/ECDF-{step00.format_filename(x)}.pdf")
            fig.savefig(figure_list[-1])
            figure_list.append(f"{directory}/ECDF-{step00.format_filename(x)}.png")
            fig.savefig(figure_list[-1])
            matplotlib.pyplot.close(fig)

        drawing_data = drawing_data.loc[(drawing_data["Score"] < 50)]

        # Hist plots
        for x in tqdm.tqdm(["Aligned length", "Maximum length"], position=1, leave=False):
            fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

            seaborn.histplot(data=drawing_data, x=x, hue="Type", palette=type_palette, stat="percent", common_norm=False, bins=100, kde=True, multiple="dodge", log_scale=True, ax=ax)

            matplotlib.pyplot.xlabel(f"{x} (bp)")
            matplotlib.pyplot.tight_layout()

            figure_list.append(f"{directory}/Hist-{step00.format_filename(x)}.pdf")
            fig.savefig(figure_list[-1])
            figure_list.append(f"{directory}/Hist-{step00.format_filename(x)}.png")
            fig.savefig(figure_list[-1])
            matplotlib.pyplot.close(fig)

        # Overlap length
        for data_type in tqdm.tqdm(list(type_palette.keys())):
            length_data[data_type] = [0 for _ in chromosome_list]

        for chromosome, data_type in tqdm.tqdm(list(itertools.product(chromosome_list, list(type_palette.keys())))):
            selected_data = drawing_data.loc[(drawing_data["Chromosome 1"] == chromosome) & (drawing_data["Type"] == data_type)]

            merged = [(0, 0)]
            for index, row in tqdm.tqdm(selected_data.iterrows(), total=len(selected_data), position=1, leave=False):
                prev_start, prev_end = merged[-1]
                start, end = row["Start 1"], row["End 1"]

                if start <= prev_end:
                    merged[-1] = (prev_start, max(prev_end, end))
                else:
                    merged.append((start, end))

            length_data.loc[chromosome, data_type] = sum(end - start for start, end in merged)
        print(length_data)

        # SD length bar
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

        matplotlib.pyplot.bar(range(len(chromosome_list)), length_data["Length"], width=0.8, align="center", color="tab:gray", label="Total")
        matplotlib.pyplot.bar(range(len(chromosome_list)), length_data["BISER"], width=-0.4, align="edge", color=type_palette["BISER"], linewidth=0, label="BISER")
        matplotlib.pyplot.bar(range(len(chromosome_list)), length_data["Reference"], width=0.4, align="edge", color=type_palette["Reference"], linewidth=0, label="Reference")

        for i, chromosome in tqdm.tqdm(enumerate(chromosome_list), total=len(chromosome_list)):
            input_length = length_data.loc[chromosome, "BISER"]
            total_length = length_data.loc[chromosome, "Length"]
            matplotlib.pyplot.text(i, input_length, f"{input_length / total_length * 100:.2f}%", color="black", fontsize="xx-small", horizontalalignment="right", verticalalignment="bottom", rotation="vertical", path_effects=step00.path_effects)

        for i, chromosome in tqdm.tqdm(enumerate(chromosome_list), total=len(chromosome_list)):
            input_length = length_data.loc[chromosome, "Reference"]
            total_length = length_data.loc[chromosome, "Length"]
            matplotlib.pyplot.text(i, input_length, f"{input_length / total_length * 100:.2f}%", color="black", fontsize="xx-small", horizontalalignment="left", verticalalignment="bottom", rotation="vertical", path_effects=step00.path_effects)

        matplotlib.pyplot.xlabel("Chromosome")
        matplotlib.pyplot.xticks(range(len(chromosome_list)), chromosome_list, fontsize="x-small", rotation="vertical")
        matplotlib.pyplot.ylabel("Length (bp)")
        matplotlib.pyplot.legend(loc="upper right")
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.tight_layout()

        figure_list.append(f"{directory}/SD-length.pdf")
        fig.savefig(figure_list[-1])
        figure_list.append(f"{directory}/SD-length.png")
        fig.savefig(figure_list[-1])
        matplotlib.pyplot.close(fig)

        # Violin plot
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

        seaborn.violinplot(data=drawing_data, x="Chromosome 1", order=chromosome_list, y="Score", hue="Type", hue_order=list(type_palette.keys()), palette=type_palette, inner="box", cut=0, legend="full", ax=ax)

        matplotlib.pyplot.xlabel("Chromosome")
        matplotlib.pyplot.xticks(fontsize="x-small", rotation="vertical")
        matplotlib.pyplot.tight_layout()

        figure_list.append(f"{directory}/Score-all.pdf")
        fig.savefig(figure_list[-1])
        figure_list.append(f"{directory}/Score-all.png")
        fig.savefig(figure_list[-1])
        matplotlib.pyplot.close(fig)

        for chromosome in tqdm.tqdm(chromosome_list):
            selected_input_data = input_data.loc[(input_data["Chromosome 1"] == chromosome)]
            selected_reference_data = reference_data.loc[(reference_data["Chromosome 1"] == chromosome)]

            input_rows = [(0, 0)]
            for index, row in tqdm.tqdm(selected_input_data.iterrows(), total=len(selected_input_data), position=1, leave=False):
                prev_start, prev_end = input_rows[-1]
                start, end = row["Start 1"], row["End 1"]

                if start <= prev_end:
                    input_rows[-1] = (prev_start, max(prev_end, end))
                else:
                    input_rows.append((start, end))

            reference_rows = [(0, 0)]
            for index, row in tqdm.tqdm(selected_reference_data.iterrows(), total=len(selected_reference_data), position=1, leave=False):
                prev_start, prev_end = reference_rows[-1]
                start, end = row["Start 1"], row["End 1"]

                if start <= prev_end:
                    reference_rows[-1] = (prev_start, max(prev_end, end))
                else:
                    reference_rows.append((start, end))

            intersections = list()
            i = j = 0
            while (i < len(input_rows)) and (j < len(reference_rows)):
                input_start, input_end = input_rows[i]
                reference_start, reference_end = reference_rows[j]

                overlap_start = max(input_start, reference_start)
                overlap_end = min(input_end, reference_end)

                if overlap_start < overlap_end:
                    intersections.append((overlap_start, overlap_end))

                if input_end <= reference_end:
                    i += 1
                else:
                    j += 1

            overlap_length = sum(end - start for start, end in intersections)

            # SD plot
            fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

            for start, end in tqdm.tqdm(intersections, position=1, leave=False):
                matplotlib.pyplot.barh(0.0, (end - start), left=start, color="tab:gray", linewidth=0)

            for index, row in tqdm.tqdm(selected_input_data.iterrows(), total=len(selected_input_data), position=1, leave=False):
                matplotlib.pyplot.barh(1.0, (row["End 1"] - row["Start 1"]), left=row["Start 1"], color=type_palette["BISER"], linewidth=0)

            for index, row in tqdm.tqdm(selected_reference_data.iterrows(), total=len(selected_reference_data), position=1, leave=False):
                matplotlib.pyplot.barh(-1.0, (row["End 1"] - row["Start 1"]), left=row["Start 1"], color=type_palette["Reference"], linewidth=0)

            matplotlib.pyplot.text(length_data.loc[chromosome, "Length",] / 2, 1.0, f"{overlap_length / length_data.loc[chromosome, 'BISER'] * 100:.2f}%", color="black", fontsize="large", horizontalalignment="center", verticalalignment="center", path_effects=step00.path_effects)
            matplotlib.pyplot.text(length_data.loc[chromosome, "Length",] / 2, 0.0, f"{overlap_length / length_data.loc[chromosome, 'Length'] * 100:.2f}%", color="black", fontsize="large", horizontalalignment="center", verticalalignment="center", path_effects=step00.path_effects)
            matplotlib.pyplot.text(length_data.loc[chromosome, "Length",] / 2, -1.0, f"{overlap_length / length_data.loc[chromosome, 'Reference'] * 100:.2f}%", color="black", fontsize="large", horizontalalignment="center", verticalalignment="center", path_effects=step00.path_effects)

            matplotlib.pyplot.title(chromosome)
            matplotlib.pyplot.xlabel("Location (bp)")
            matplotlib.pyplot.ylabel("Data type")
            matplotlib.pyplot.yticks([-1.0, 0.0, 1.0], ["Reference", "Overlap", "BISER"])
            matplotlib.pyplot.grid(True)
            matplotlib.pyplot.tight_layout()

            figure_list.append(f"{directory}/Overlap-{chromosome}.pdf")
            fig.savefig(figure_list[-1])
            figure_list.append(f"{directory}/Overlap-{chromosome}.png")
            fig.savefig(figure_list[-1])
            matplotlib.pyplot.close(fig)

        for query_chromosome, target_chromosome in tqdm.tqdm(sorted(chromosome_pair_set)):
            drawing_data = pandas.concat([input_data.loc[(input_data["Chromosome 1"] == query_chromosome) & (input_data["Chromosome 2"] == target_chromosome)], reference_data.loc[(reference_data["Chromosome 1"] == query_chromosome) & (reference_data["Chromosome 2"] == target_chromosome)]], join="inner", ignore_index=True)
            drawing_data = drawing_data.loc[(drawing_data["Score"] < 50)]

            # Scatter
            fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

            seaborn.scatterplot(data=drawing_data, x="Start 1", y="Start 2", hue="Score", palette="YlOrRd_r", style="Type", markers=type_marker, legend="brief", s=100, edgecolor=None, rasterized=True, ax=ax)

            matplotlib.pyplot.xlabel(f"Query: {query_chromosome} (bp)")
            matplotlib.pyplot.ylabel(f"Target: {target_chromosome} (bp)")
            matplotlib.pyplot.title(f"{query_chromosome}→{target_chromosome}")
            matplotlib.pyplot.legend(loc="lower right")
            matplotlib.pyplot.tight_layout()

            figure_list.append(f"{directory}/Scatter-{query_chromosome}-{target_chromosome}.pdf")
            fig.savefig(figure_list[-1])
            figure_list.append(f"{directory}/Scatter-{query_chromosome}-{target_chromosome}.png")
            fig.savefig(figure_list[-1])
            matplotlib.pyplot.close(fig)

        with zipfile.ZipFile(args.output, mode="w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zip_file:
            for figure in tqdm.tqdm(figure_list):
                zip_file.write(figure, arcname=os.path.basename(figure))
