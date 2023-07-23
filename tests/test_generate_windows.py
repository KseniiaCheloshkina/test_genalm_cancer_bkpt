import sys
import os

sys.path.append(os.getcwd())
from src.generate_windows import get_sequence


def test_get_sequence():
    assert get_sequence(chrom="1", start="40000", end="40010") == "gcctcatgga"
