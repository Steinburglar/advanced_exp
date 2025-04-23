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


def load_AB_txt(rel_path):
    """
    loads AB coincidence data from txt file
    """
    return

def load_muon_decay_txt(rel_path):
    """
    loads muon decay data from txt file
    """
    df = pd.read_csv(rel_path, sep="\t", header=None, names=["time (ns)", "count"])
    df = df.set_index("time (ns)")
    df = df.iloc[10:]
    
    return df