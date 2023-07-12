import json
import numpy as np


def write_json(fname, data):
    with open(fname, 'w') as f:
        json.dump(data, f)
        

def read_json(fname):
    with open(fname, 'r') as f:
        recv = json.load(f)
    return recv


def get_begin_pos(starts):
    eq = list(map(lambda i: starts[i] != starts[i - 1] if i > 0 else True, range(len(starts))))
    begin_pos = np.where(eq)[0].astype(int)
    return begin_pos


def get_slice(x, begin_pos):
    slices = list(map(
        lambda i: x[begin_pos[i] : begin_pos[i + 1]] if i < len(begin_pos) - 1 else x[begin_pos[i] : len(x)],
        range(len(begin_pos))
    ))
    return slices
