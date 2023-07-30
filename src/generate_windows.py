""" Generate window arounf breakpoint or random location (positive and negative examples respectively)"""
import argparse
import json
import os
from typing import Tuple
import requests
import pandas as pd
import tqdm
import numpy as np


def generate_negative_different_length(win_len: int, negative_path: str) -> pd.DataFrame:
    df = pd.read_csv(negative_path)[['position', 'chromosome', 'label']]
    all_seq = []
    for _, bkpt in tqdm.tqdm(df.iterrows()):
        win_start, win_end = generate_window(
            chrom=bkpt["chromosome"], pos=int(bkpt["position"]), win_len=win_len
        )
        all_seq.append((win_start, win_end))
    df["win_start"] = [el[0] for el in all_seq]
    df["win_end"] = [el[1] for el in all_seq]   
    return df
    
    
def generate_window(chrom: str, pos: int, win_len: int) -> Tuple[int, int]:
    """ Generates window around point taking into account chromosome lengths

    Args:
        chrom (str): chromosome number
        pos (int): position (coordinate)
        win_len (int): window length to generate (total nucleotides)

    Returns:
        Tuple[int, int]: start-end of window
    """
    # to check if chromosome is shorter
    with open("data/chr_lengths.json", "r") as f:
        chr_lengths = json.load(f)
    start = max(0, pos - round(win_len / 2))
    end = min(pos + round(win_len / 2) - 1, chr_lengths[str(chrom)])
    return start, end


def get_sequence(chrom: str, start: str, end: str) -> str:
    """ Get DNA sequence by coordinates using API. Good solution only for small number of requests
    It takes a lot of time on a big dataset...
    Args:
        chrom (str): chromosome number
        start (str): start position of the window
        end (str): end position of the window

    Returns:
        str: DNA sequence
    """
    addr = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{chrom};start={start};end={end}"
    answ = requests.post(addr).json()["dna"]
    return answ


def get_positive_windows(csv_path: str, win_len: int) -> pd.DataFrame:
    """ For each breakpoint in a file generate window around it and set positive label

    Args:
        csv_path (str): path to breakpoints data
        win_len (int): window length

    Returns:
        pd.DataFrame: resulting DF
    """
    df = pd.read_csv(csv_path)
    all_seq = []
    for _, bkpt in tqdm.tqdm(df.iterrows()):
        win_start, win_end = generate_window(
            chrom=bkpt["chr"], pos=int(bkpt["start"]), win_len=win_len
        )
        all_seq.append((win_start, win_end))
    df["win_start"] = [el[0] for el in all_seq]
    df["win_end"] = [el[1] for el in all_seq]   
    df = df[["chr", "start", "win_start", "win_end", "cancer_type"]]
    df = df.rename(columns={"start": "position", "chr": "chromosome"})
    df["label"] = 1
    return df


def get_negative_windows(n_points: int, win_len: int) -> pd.DataFrame:
    """ Generates specified number of negative examples randomly from each chromosome.
    For each chromosome points are generated uniformly based on its length and excluding telomeres

    Args:
        n_points (int): number of points to generate
        win_len (int): window length

    Returns:
        pd.DataFrame: resulting dataframe
        ,
        
    """
    with open("data/chr_lengths.json", "r") as f:
        chr_lengths = json.load(f)
    mean_telomeres_len = 10000
    num_points_per_chr = round(n_points / len(chr_lengths))
    all_points = []
    for chrom, chr_length in tqdm.tqdm(chr_lengths.items()):
        positions = np.random.uniform(
            low=mean_telomeres_len,
            high=chr_length - mean_telomeres_len,
            size=num_points_per_chr,
        ).astype(int)
        df = pd.DataFrame(data=positions, columns=["position"])
        df["chromosome"] = chrom
        values = []
        for pos in positions:
            win_start, win_end = generate_window(chrom, pos, win_len)
            values.append((win_start, win_end))
        df["win_start"] = [el[0] for el in values]
        df["win_end"] = [el[1] for el in values]
        all_points.append(df)
    df_all = pd.concat(all_points)
    df_all = df_all.sort_values(['chromosome', 'win_start'])
    df_all["label"] = 0
    return df_all



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--win_len", help="length of window to generate", default=512, type=int
    )
    parser.add_argument(
        "--run_number", help="order number of window length", default=1, type=int
    )
    args = parser.parse_args()
    main_path = "data/dataset/"
    # save positive
    df_pos = get_positive_windows(csv_path="data/breakpoints_wo_bad_regions.csv", 
                                  win_len=args.win_len)
    df_pos.to_csv(f"{main_path}positive_all_cancers_{args.win_len}.csv")
    
    # ATTENTION: do it only 1 time for 1 window length. Datasets with all the rest window lengths
    # should contains the same set of negative points (and not generating a differet set)
    # generate_negative - for the first time
    neg_path = f"{main_path}negative_all_cancers_{args.win_len}.csv"
    if args.run_number == 1:
        print("generate new negatives")
        df_neg = get_negative_windows(n_points=1000000, win_len=args.win_len)
    else:
        print("expand existing negatives")
        # use existing negative set and expand window length
        df_neg = generate_negative_different_length(
            win_len=args.win_len, 
            negative_path=f"{main_path}negative_all_cancers_512.csv")
    df_neg.to_csv(neg_path)
    
    # save to bed format to finally get DNA sequence
    flnms = [fl for fl in os.listdir(main_path) if fl.endswith(".csv")]
    for fl in flnms:
        df_windows = pd.read_csv(os.path.join(main_path, fl))
        df_windows['chromosome'] = df_windows['chromosome'].map(lambda x: 'chr' + str(x))
        df_windows[['chromosome', 'win_start', 'win_end']].to_csv(
            os.path.join(main_path, fl.replace(".csv", ".bed")), 
            sep="\t",
            header=False,
            index=False
        )
