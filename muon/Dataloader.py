#dataloader for muon decay experiment

import pandas as pd
import math
import numpy as np
from cavendish.utils.Functions import *
from millikan.functions import *
from compton.Functions import *
from compton.Analysis import *
from compton.Dataloader import *
from positron.Functions import *
from positron.Analysis import *
from positron.Dataloader import *


def load_ABT_txt(rel_path):
    """
    loads AB coincidence data from txt file
    """
    df = pd.read_csv(rel_path, sep="\t")
    df.set_index("Voltage", inplace=True)
    return df

def load_muon_decay_txt(rel_path):
    """
    loads muon decay data from txt file
    """
    df = pd.read_csv(rel_path, sep="\t", header=None, names=["Channel", "count"])
    df["time (microseconds)"] = df["Channel"] * (8/1000)
    df = df.set_index("time (microseconds)")
    df.drop(columns=["Channel"], inplace=True)
    df = df.iloc[10:]
    
    return df