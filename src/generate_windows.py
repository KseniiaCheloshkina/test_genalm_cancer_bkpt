import argparse
import json


def generate_window(chr: str, pos: int, win_len: int):
    # to check if chromosome is shorter
    with open("data/chr_lengths.json", "r") as f:
        chr_lengths = json.load(f)
    start = max(0, pos - round(win_len/2))
    end = min(pos + round(win_len/2) - 1, chr_lengths[chr])
    return start, end


def get_sequence():
    pass


def get_positive_windows():
    # output: cancer type, chromosome, win_start, win_end, breakpoint position, win_sequence, label
    pass


def get_negative_windows():
    # output: cancer type, chromosome, win_start, win_end, breakpoint position, win_sequence, label
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--win_len", help="length of window to generate", default=512, type=int)
    args = parser.parse_args()
    # generate_window(win_len=args.win_len)
