import argparse
import json
import os
from typing import Tuple
import requests
import pandas as pd
import tqdm
import numpy as np


def generate_window(chrom: str, pos: int, win_len: int) -> Tuple[int, int]:
    # to check if chromosome is shorter
    with open("data/chr_lengths.json", "r") as f:
        chr_lengths = json.load(f)
    start = max(0, pos - round(win_len / 2))
    end = min(pos + round(win_len / 2) - 1, chr_lengths[str(chrom)])
    return start, end


def get_sequence(chrom: str, start: str, end: str) -> str:
    # It takes a lot of time ona big dataset...
    addr = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{chrom};start={start};end={end}"
    answ = requests.post(addr).json()["dna"]
    return answ


def get_positive_windows(csv_path: str, win_len: int) -> pd.DataFrame:
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
    args = parser.parse_args()
    
    main_path = "data/dataset/"
    # save positive
    df_pos = get_positive_windows(csv_path="data/breakpoints_wo_bad_regions.csv", 
                                  win_len=args.win_len)
    df_pos.to_csv(f"{main_path}positive_all_cancers_{args.win_len}.csv")
    # generate_negative
    df_neg = get_negative_windows(n_points=1000000, win_len=args.win_len)
    df_neg.to_csv(f"{main_path}negative_all_cancers_{args.win_len}.csv")
    
    # save to ved format to finally get DNA sequence
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
