import os
import pandas as pd 
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

from mixcurve import MixingConverter

class ProbabilisticMixingConverter(MixingConverter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
