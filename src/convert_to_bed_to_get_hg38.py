""" Script to prepare breakpoints dataset for conversion and process the results"""
import pandas as pd
from typing import List


def create_bed_format(df_row: pd.DataFrame) -> str:
    """Convert one row to bed format

    Args:
        df_row (pd.DataFrame): dataframe with breakpoints data

    Returns:
        str: breakpoint location in bed format
    """
    if df_row["chr_range"] <= 10:
        return (
            "chr"
            + df_row["chr"]
            + ":"
            + str(df_row["chr_bkpt"])
            + "-"
            + str(df_row["chr_bkpt"] + 1)
        )
    else:
        return "bad"


def process_initial_data(filepath: str) -> pd.DataFrame:
    """Load initial breakpoints data and remove some regions:
        * X, Y chromosome
        * breakpoints with low precision

    Args:
        filepath (str): path to the file

    Returns:
        pd.DataFrame: resulting breakpoints
    """
    all_bkpt = pd.read_csv(filepath)
    print("Total:", all_bkpt.shape[0])  # 487 425
    print(all_bkpt["chr"].unique())
    all_bkpt = all_bkpt[~all_bkpt["chr"].isin(["X", "Y"])]
    all_bkpt["chr_bkpt"] = all_bkpt["chr_bkpt"].astype(int)

    all_bkpt["formatted"] = all_bkpt.apply(create_bed_format, axis=1)
    good_bkpt = all_bkpt[all_bkpt["formatted"] != "bad"]
    print("After filtration:", good_bkpt.shape[0])  # 472 573
    return good_bkpt


def prepare_file_for_convertation(filepath: str, out_file: str) -> None:
    """Process initial breakpoints and save in txt format

    Args:
        filepath (str): path to the initial file
        out_file (str): path to the output txt file
    """
    good_bkpt = process_initial_data(filepath)

    with open(out_file, "w") as f:
        for _, row in good_bkpt.iterrows():
            f.write(row["formatted"] + "\n")


def get_final_file(init_file: str, err_file: str, result_file: str) -> List[str]:
    """Get correct hg38 coordinates taking into account errors list

    Args:
        init_file (str): list of hg19 coordinates from function `prepare_file_for_convertation`
        err_file (str): file with conversion errors
        result_file (str): result file with hg38 coordinates

    Returns:
        List[str]: hg38 coodinates with inserted errors
    """
    hg19_data = []
    with open(init_file, "r") as f:
        for line_hg19 in f.readlines():
            hg19_data.append(line_hg19.replace("\n", ""))

    # get all hg19 errors
    all_err = []
    with open(err_file, "r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                all_err.append(line.replace("\n", ""))

    all_err_indices = []
    err_order_num = 0
    for ind_hg19, line_hg19 in enumerate(hg19_data):
        if line_hg19 == all_err[err_order_num]:
            all_err_indices.append(ind_hg19)
            if len(all_err) > err_order_num + 1:
                err_order_num += 1

    # get all hg38 points
    hg38_data = []
    with open(result_file, "r") as f:
        for line in f.readlines():
            hg38_data.append(line.replace("\n", ""))

    for el in all_err_indices:
        hg38_data.insert(el, "bad coordinate")

    return hg38_data


def extract_chr(x: str) -> str:
    return (
        x.split(":")[0].replace("chr", "")
        if x != "bad coordinate"
        else "bad coordinate"
    )


def extract_position(x: str) -> str:
    return x.split(":")[1].split("-")[0] if x != "bad coordinate" else "bad coordinate"


def get_full_dataset(
    initial_bkpt_path: str,
    init_file: str,
    err_file: str,
    result_file: str,
    out_file: str,
) -> None:
    """Add as columns coordinates of hg38 in source dataset

    Args:
        initial_bkpt_path (str): initial file with breakpoints data
        init_file (str): list of hg19 coordinates from function `prepare_file_for_convertation`
        err_file (str): file with conversion errors
        result_file (str): result file with hg38 coordinates
        out_file (str): output file
    """
    good_bkpt = process_initial_data(initial_bkpt_path)
    hg38_coord = get_final_file(init_file, err_file, result_file)
    hg38_coord = pd.DataFrame(hg38_coord, columns=["coordinate"])
    hg38_coord["hg38_chr"] = hg38_coord["coordinate"].map(extract_chr)
    hg38_coord["hg38_coord"] = hg38_coord["coordinate"].map(extract_position)
    assert hg38_coord.shape[0] == good_bkpt.shape[0]
    good_bkpt = pd.concat(
        [good_bkpt.reset_index(drop=True), hg38_coord.reset_index(drop=True)], axis=1
    )
    print(good_bkpt[["chr", "chr_bkpt", "hg38_chr", "hg38_coord"]].head())
    print((good_bkpt["hg38_chr"] != good_bkpt["chr"]).astype(int).sum())
    good_bkpt[good_bkpt["hg38_chr"] != good_bkpt["chr"]].to_csv("data/bad_rows.csv")
    good_bkpt = good_bkpt[good_bkpt["hg38_chr"] == good_bkpt["chr"]]
    good_bkpt.to_csv(out_file)


if __name__ == "__main__":
    initial_bkpt_path = "../cancer_breakpoints_hotspots_prediction/data/raw breakpoints/all_cancer_data_eda.csv"

    ## STEP 1.
    # create file to convert to hg38
    # prepare_file_for_convertation(
    #     filepath=initial_bkpt_path,
    #     out_file="data/hg19_breakpoints.txt"
    # )

    # then upload this file to site https://genome.ucsc.edu/cgi-bin/hgLiftOver
    # Result: Successfully converted 472429 records. Conversion failed on 144 records

    ## STEP 2.
    # merge error file
    get_full_dataset(
        initial_bkpt_path=initial_bkpt_path,
        init_file="data/hg19_breakpoints.txt",
        err_file="data/hg19_breakpoints_err.txt",
        result_file="data/hg38_breakpoints.bed",
        out_file="data/hg38_breakpoints_wo_err.csv",
    )
