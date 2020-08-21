from __future__ import division
import os
from collections import Counter

def get_path():
    return os.path.expanduser("~/GitHub/Phylo_Evol_Timeseries")

def get_to_keep():
    return ['SNP', 'INS', 'DEL']

def common_entries(*dcts):
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)

def get_bacillus_mut_bias():
    return  1.2739

def get_bacillus_mut_rate():
    return 3.35 * (10**-10)

def get_bacillus_indel_rate():
    return 1.20 * (10**-10)

def get_colors():
    return {'0':'#87CEEB', '1': '#FFA500', '2':'#FF6347'}


def get_ts_tv_dict():

    ts_tv_dict = {
    ("A", "C"):"V", ("A", "G"):"S", ("A", "T"):"V",
    ("C", "A"):"V", ("C", "G"):"V", ("C", "T"):"S",
    ("G", "A"):"S", ("G", "C"):"V", ("G", "T"):"V",
    ("T", "A"):"V", ("T", "C"):"S", ("T", "G"):"V"}

    return ts_tv_dict
