import pandas as pd 
import numpy as np
from scipy.optimize import curve_fit

from mixcurve import MixingConverter

class ModelException(Exception):
    pass

class NotImplementedException(Exception):
    pass

class ProbabilisticMixingConverter(MixingConverter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_loaded = False


    def _fwd_func(self, x, a, c, d):
        return a*np.exp(-c*x)+d

    def _percent_corrected(self, straight, mix, pool, target):
        try:
            return 100. * (mix - straight) / (target - straight)
        except ZeroDivisionError:
            return 100

    def build_model(self, sigma_seconds=0.5, sigma_percent=5, n_iters=1000):
        if self.method != "curve_fit":
            raise ModelException("Model building is supported for the "
                                 "`curve_fit` MixingConverter only!")
        s = len(self.table)
        t = self.table
        out = np.zeros([n_iters,3])
        for i in range(n_iters):
            x = self.table.PTT + np.random.normal(0, sigma_seconds, size=s)
            y = self.table.factor_level + np.random.normal(0, sigma_percent, size=s)
            params, _ = curve_fit(self._fwd_func, x, y,p0=(40000, 0.1, 20), ftol=0.001, xtol=0.001)
            out[i, :] = params

        self.out = out
        self.sigma_seconds = sigma_seconds
        self.sigma_percent = sigma_percent
        self.model_loaded = True

    def save_model(self, filename):
        raise NotImplementedException()

    def load_model(self, filename):
        raise NotImplementedException()

    def pt_percent_corrected(self, pt_straight, pt_mix, pt_pool, n_iters=10e4):
        if not self.model_loaded:
            raise ModelException("Build or load a model first")
        raise NotImplementedException()
        
    def ptt_percent_corrected(self, ptt_straight, ptt_mix, ptt_pool, n_iters=100000):
        if not self.model_loaded:
            raise ModelException("Build or load a model first")

        # Step 1: Draw random values from normal distribution for PTT/Mix/Pool
        # The distribution represents the error inherent to the measured results
        assay_params = {"sigma": self.sigma_seconds, "n_iters": n_iters}
        PTT = self._normal_dist(ptt_straight, **assay_params)
        Mix = self._normal_dist(ptt_mix, **assay_params)
        Pool = self._normal_dist(ptt_pool, **assay_params)

        # Step 2: Select random sets of parameters from the curve_fit model
        curve_fit_param_indices = np.random.randint(0, self.out.shape[0]-1, size=n_iters)
        curve_fit_params = self.out[curve_fit_param_indices, :].T

        # Step 3: Transform the PTT/Mix/Pool values (seconds) to factor level, 
        # based on the curve_fit_params selected in step 2
        # _fl == "_factor_level"
        PTT_fl = self._fwd_func(PTT, *curve_fit_params)
        Mix_fl = self._fwd_func(Mix, *curve_fit_params)
        Pool_fl = self._fwd_func(Pool, *curve_fit_params)

        # Step 4: Calculate the target factor level (Target_fl)
        # And then the perecent_corrected
        Target_fl = (PTT_fl + Pool_fl) / 2.0
        percent_corrected = self._percent_corrected(PTT_fl, Mix_fl, Pool_fl, Target_fl)

        return percent_corrected

    def _normal_dist(self, value, sigma, n_iters):
        return np.random.normal(loc=value, scale=sigma, size=n_iters)