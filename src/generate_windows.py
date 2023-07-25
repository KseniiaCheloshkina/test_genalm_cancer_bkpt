import argparse
import json
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
    addr = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{chrom};start={start};end={end}"
    answ = requests.post(addr).json()["dna"]
    # TODO: check if its needed to set it to lower or upper case
    return answ


def get_positive_windows(csv_path: str, win_len: int) -> None:
    df = pd.read_csv(csv_path)
    all_seq = []
    for _, bkpt in tqdm.tqdm(df.iterrows()):
        win_start, win_end = generate_window(
            chrom=bkpt["chr"], pos=int(bkpt["start"]), win_len=win_len
        )
        # TODO: fix
        # dna_seq = get_sequence(
        #     chrom=bkpt["chr"], start=str(win_start), end=str(win_end)
        # )
        dna_seq = "a"
        all_seq.append((win_start, win_end, dna_seq))
    df["win_start"] = [el[0] for el in all_seq]
    df["win_end"] = [el[1] for el in all_seq]
    df["dna_seq"] = [el[2] for el in all_seq]    
    df.to_csv(csv_path.replace(".csv", f"_{win_len}.csv"))
    df = df[["chr", "start", "win_start", "win_end", "dna_seq", "cancer_type"]]
    df = df.rename(columns={"start": "position", "chr": "chromosome"})
    df["label"] = 1
    print("With N: ", df[df["sequence"].str.contains("N")].shape[0])    
    return df


def get_negative_windows(n_points: int, win_len: int) -> pd.DataFrame:
    # output: chromosome, win_start, win_end, position, dna_seq, label
    with open("data/chr_lengths.json", "r") as f:
        chr_lengths = json.load(f)
    mean_telomeres_len = 10000
    num_points_per_chr = round(n_points / len(chr_lengths))
    all_points = []
    for chrom, chr_length in chr_lengths.items():
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
            dna_seq = get_sequence(chrom=chrom, start=str(win_start), end=str(win_end))
            values.append((win_start, win_end, dna_seq))
        df["win_start"] = [el[0] for el in values]
        df["win_end"] = [el[1] for el in values]
        df["dna_seq"] = [el[2] for el in values]
        all_points.append(df)
    df_all = pd.concat(all_points)
    df_all["label"] = 0
    return df_all


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--win_len", help="length of window to generate", default=512, type=int
    )
    args = parser.parse_args()
    df_new = get_positive_windows(csv_path="data/breakpoints_wo_bad_regions.csv", win_len=args.win_len)
    # df_new = get_negative_windows(n_points=40, win_len=args.win_len)
    # print(df_new.head())
    # df_new.to_csv("data/tmp.csv")
