""" Tests for file src/filter_bad_breakpoints.py"""
import sys
import os
import pandas as pd

sys.path.append(os.getcwd())
from src.filter_bad_breakpoints import get_intersected_rows


def test_get_intersected_rows():
    test_df_1 = [
        {'chr': '1', 'start': 0, 'end': 10000},
        {'chr': '1', 'start': 10000, 'end': 10005},
        {'chr': '2', 'start': 20000, 'end': 30000},
        {'chr': '2', 'start': 30000, 'end': 40000}
        ]
    test_df_2 = [
        {'chr': '1', 'start': 500, 'end': 10500},
        {'chr': '1', 'start': 0, 'end': 10005},
        {'chr': '2', 'start': 15000, 'end': 31000},
        {'chr': '2', 'start': 25000, 'end': 28000}
        ]
    df_intersected = get_intersected_rows(pd.DataFrame(test_df_1), pd.DataFrame(test_df_2))
    assert df_intersected.shape[0] == 7
    assert df_intersected['start_1'].sum() == 56000
