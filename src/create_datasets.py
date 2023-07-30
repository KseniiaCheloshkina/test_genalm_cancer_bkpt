""" Script to collect meta and sequence data into one file,
merge positive and negative examples, and split by cancer type"""
import os
import sys
import argparse
import re
import tqdm
import pandas as pd
from sklearn.utils import shuffle

sys.path.append(os.getcwd())
from src.filter_bad_breakpoints import get_intersected_rows
from src.generate_windows import generate_window


def merge_meta_and_seq(meta_path: str, seq_path: str) -> pd.DataFrame:
    """ Collects meta data and sequence data into one dataframe

    Args:
        meta_path (str): path to data with meta information
        seq_path (str): path to sequence data

    Returns:
        pd.DataFrame: resulting dataframe
    """
    df_meta = pd.read_csv(meta_path, dtype=str).drop(["Unnamed: 0"], axis=1)
    all_seq = []
    with open(seq_path, "r") as f:
        for line in f.readlines():
            position, dna_seq = line.split("\t")
            cur_point = {
                "chr": position.split(":")[0].replace("chr", ""),
                "start": position.split(":")[1].split("-")[0],
                "end": position.split(":")[1].split("-")[1],
                "dna_seq": re.sub("\n$", "", dna_seq),
            }
            all_seq.append(cur_point)
    df_meta_seq = pd.DataFrame(all_seq)
    df_meta_all = pd.concat([df_meta, df_meta_seq], axis=1)
    assert df_meta_all[df_meta_all["chromosome"] != df_meta_all["chr"]].shape[0] == 0
    assert df_meta_all[df_meta_all["win_start"] != df_meta_all["start"]].shape[0] == 0
    assert df_meta_all[df_meta_all["win_end"] != df_meta_all["end"]].shape[0] == 0
    return df_meta_all


def prepare_data(
    pos_path: str,
    pos_path_seq: str,
    neg_path: str,
    neg_path_seq: str
):
    # read positive and merge sequence and meta data
    df_pos = merge_meta_and_seq(meta_path=pos_path, seq_path=pos_path_seq)
    df_pos = df_pos[
        ["chr", "start", "end", "position", "cancer_type", "dna_seq", "label"]
    ]
    # read negative and merge sequence and meta data
    df_neg = merge_meta_and_seq(meta_path=neg_path, seq_path=neg_path_seq)
    df_neg = df_neg[["chr", "start", "end", "position", "dna_seq", "label"]]
    # remove bad regions from negatives
    df_bad_regions = pd.read_csv("data/all_excluded_regions.csv").rename(
        columns={"chrom": "chr", "chromStart": "start", "chromEnd": "end"}
    )
    df_bad_regions['chr'] = df_bad_regions['chr'].map(lambda x: x.replace("chr", ""))
    df_bad_regions['chr'] = df_bad_regions['chr'].astype(str)
    df_intersected = get_intersected_rows(df_neg, df_bad_regions)
    df_neg = pd.merge(df_neg, df_intersected, how="left")
    df_neg = df_neg[df_neg["chr_1"].isnull()]
    df_neg.drop(['chr_1', 'start_1', 'end_1', 'Unnamed: 0'], axis=1, inplace=True)
    df_neg = df_neg.drop_duplicates()
    print("Negative examples after removal of excluded regions", df_neg.shape[0])
    return df_pos, df_neg


def get_dataset_for_cancer_type(
    pos_path: str,
    pos_path_seq: str,
    neg_path: str,
    neg_path_seq: str,
    out_folder: str,
    n_times_neg_more: int,
    win_len: int
) -> None:
    """ Saves final dataset for training a model:
    * prepares one file per cancer type
    * joins positive and negative examples
    * removes their intersections
    * transform DNA sequence to upper case (as models vocabularies do not contain lowercase )

    Args:
        pos_path (str): path to meta data for positive examples
        pos_path_seq (str): path to sequence data for positive examples
        neg_path (str): path to meta data for negative examples
        neg_path_seq (str): path to sequence data for negative examples
        out_folder (str): folder to save results
        n_times_neg_more (int): The class balance in each dataset will be 
            1:`n_times_neg_more` (positive: negative).
        win_len (int): Window length (used to name file)
    """
    df_pos, df_neg = prepare_data(pos_path=pos_path, pos_path_seq=pos_path_seq, neg_path=neg_path, neg_path_seq=neg_path_seq)
    # split by cancer type
    cancers = df_pos["cancer_type"].unique()
    for cancer_type in tqdm.tqdm(cancers):
        df_pos_cancer = df_pos[df_pos["cancer_type"] == cancer_type].drop(
            ["cancer_type"], axis=1
        )
        # remove intersecting windows
        df_intersected = get_intersected_rows(
            df_neg[["chr", "start", "end"]], df_pos_cancer[["chr", "start", "end"]]
        )
        df_intersected = df_intersected[["chr", "start", "end"]].drop_duplicates()
        df_intersected["intersected"] = 1
        df_neg_merged = pd.merge(
            df_neg, df_intersected, on=["chr", "start", "end"], how="left"
        )
        df_neg_merged = df_neg_merged[df_neg_merged["intersected"].isnull()].drop(['intersected'], axis=1)
        # sample same number of negatives (equally distributed by chromosomes)
        n_points_per_chr = n_times_neg_more * (round(df_pos_cancer.shape[0] / 23) + 1)
        all_neg_for_cancer = []
        for chr_num in df_neg_merged["chr"].unique():
            df_neg_merged_cur = df_neg_merged[df_neg_merged["chr"] == chr_num]
            indices = df_neg_merged_cur.index.tolist()
            ind_to_take = indices[
                0 :: round(df_neg_merged_cur.shape[0] / n_points_per_chr)
            ]
            all_neg_for_cancer.append(df_neg_merged_cur.loc[ind_to_take])
        df_neg_for_cancer = pd.concat(all_neg_for_cancer)
        df_final = pd.concat([df_pos_cancer, df_neg_for_cancer], axis=0)
        # to upper case
        df_final["dna_seq"] = df_final["dna_seq"].map(lambda x: x.upper())
        df_final = shuffle(df_final)
        # save csv
        df_final.reset_index(drop=True).to_csv(
            os.path.join(out_folder, f"{cancer_type}_{n_times_neg_more}_{win_len}.csv"),
            index=False
        )


def match_similar_negatives(
    pos_path: str,
    pos_path_seq: str,
    neg_path: str,
    neg_path_seq: str,
    out_folder: str,
    n_times_neg_more: int,
    win_len: int
):
    df_pos, df_neg = prepare_data(
        pos_path=pos_path,
        pos_path_seq=pos_path_seq,
        neg_path=neg_path,
        neg_path_seq=neg_path_seq
    )
    # split by cancer type
    cancers = df_pos["cancer_type"].unique()
    for cancer_type in tqdm.tqdm(cancers):
        print(cancer_type)
        df_pos_cancer = df_pos[df_pos["cancer_type"] == cancer_type].drop(
            ["cancer_type"], axis=1
        )
        # read previous negatives for this cancer types
        df_neg_old = pd.read_csv(f"data/dataset/final/{cancer_type}_{n_times_neg_more}_512.csv")
        df_neg_old = df_neg_old[df_neg_old['label'] == 0][['chr', 'position']]
        # generate left window boundary
        all_starts = []
        for _, row in df_neg_old.iterrows():
            all_starts.append(generate_window(row['chr'], row['position'], win_len)[0])
        df_neg_old['start'] = all_starts
        df_neg_old['start'] = df_neg_old['start'].astype(str)
        # merge with current negatives
        df_neg_for_cancer = pd.merge(df_neg, df_neg_old, on=['chr', 'start'], how="inner")
        print(df_neg_for_cancer.shape[0])
        print(df_pos_cancer.shape[0])
        df_final = pd.concat([df_pos_cancer, df_neg_for_cancer], axis=0)
        # to upper case
        df_final["dna_seq"] = df_final["dna_seq"].map(lambda x: x.upper())
        df_final = shuffle(df_final)
        # save csv
        df_final.reset_index(drop=True).to_csv(
            os.path.join(out_folder, f"{cancer_type}_{n_times_neg_more}_{win_len}.csv"),
            index=False
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n_times_neg_more", help="""
        There will be `n_time` more negative examples in a dataset for each cancer types than positives
        """, default=1, type=int
    )
    parser.add_argument(
        "--win_len", help="length of window", default=512, type=int
    )
    parser.add_argument(
        "--run_number", help="order number of window length", default=1, type=int
    )
    args = parser.parse_args()
    main_input_path = "data/dataset/"
    if args.run_number == 1:
        print('generate new')
        get_dataset_for_cancer_type(
            pos_path=f"{main_input_path}positive_all_cancers_{args.win_len}.csv",
            pos_path_seq=f"{main_input_path}pos_{args.win_len}.bed",
            neg_path=f"{main_input_path}negative_all_cancers_{args.win_len}.csv",
            neg_path_seq=f"{main_input_path}neg_{args.win_len}.bed",
            out_folder=main_input_path + "final",
            n_times_neg_more=args.n_times_neg_more,
            win_len=args.win_len
        )
    else:
        print("use negatives from 512 window length")
        match_similar_negatives(
            pos_path=f"{main_input_path}positive_all_cancers_{args.win_len}.csv",
            pos_path_seq=f"{main_input_path}pos_{args.win_len}.bed",
            neg_path=f"{main_input_path}negative_all_cancers_{args.win_len}.csv",
            neg_path_seq=f"{main_input_path}neg_{args.win_len}.bed",
            out_folder=main_input_path + "final",
            n_times_neg_more=args.n_times_neg_more,
            win_len=args.win_len
        )
