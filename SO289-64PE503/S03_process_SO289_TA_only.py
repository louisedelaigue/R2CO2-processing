import numpy as np, pandas as pd, koolstof as ks, calkulate as calk
import copy, itertools
from pandas.tseries.offsets import DateOffset
import matplotlib.dates as mdates
from matplotlib import pyplot as plt

# Import logfile and dbs file
logfile = ks.read_logfile(
    "data/TA_ONLY/logfile.bak",
    methods=[
        "3C standard AT only",
    ],
)

dbs = ks.read_dbs("data/TA_ONLY/SO289_TA_only.dbs", logfile=logfile)