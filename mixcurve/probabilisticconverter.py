import pandas as pd 
import numpy as np
import pickle
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

    def _fwd_func_linear(self, x):
        return np.interp(x, self.table[self.measure], self.table.factor_level)

    def _percent_corrected(self, straight, mix, pool, target):
        try:
            return 100. * (mix - straight) / (target - straight)
        except ZeroDivisionError:
            return 100

    def build_model(self, measure, sigma_seconds=0.5, sigma_factor_level=5, n_iters=1000):
        if self.method != "curve_fit":
            raise ModelException("Model building is supported for the "
                                 "`curve_fit` MixingConverter only!")
        if measure.upper() not in ["PTT", "PT"]:
            raise ModelException("`measure` must be 'PT' or 'PTT'")
        s = len(self.table)
        t = self.table
        out = np.zeros([n_iters,3])
        for i in range(n_iters):
            x = self.table[measure] + np.random.normal(0, sigma_seconds, size=s)
            y = self.table.factor_level + np.random.normal(0, sigma_factor_level, size=s)
            params, _ = curve_fit(self._fwd_func, x, y,p0=(40000, 0.1, 20), ftol=0.001, xtol=0.001)
            out[i, :] = params

        self.out = out
        self.measure = measure
        self.sigma_seconds = sigma_seconds
        self.sigma_factor_level = sigma_factor_level
        self.model_loaded = True

    def save_model(self, filename):
        model = {
            "out": self.out,
            "sigma_seconds": self.sigma_seconds,
            "sigma_factor_level": self.sigma_factor_level,
            "model_loaded": self.model_loaded,
            "measure": self.measure
        }
        with open(filename, 'wb') as f:
            pickle.dump(model, f)

    def load_model(self, filename):
        with open(filename, 'rb') as f:
            model = pickle.load(f)
        self.out = model["out"]
        self.sigma_seconds = model["sigma_seconds"]
        #self.sigma_factor_level = model["sigma_factor_level"]
        self.model_loaded = model["model_loaded"]
        self.measure = model["measure"]

        
    def percent_corrected_dist(self, straight, mix, pool, n_iters=100000, sigma=None, use_method='curve'):
        if not self.model_loaded:
            raise ModelException("Build or load a model first")

        # Step 1: Draw random values from normal distribution for Stright/Mix/Pool
        # The distribution represents the error inherent to the measured results
        if sigma:
            self.sigma_seconds = sigma
        assay_params = {"sigma": self.sigma_seconds, "n_iters": n_iters}
        straight_dist = self._normal_dist(straight, **assay_params)
        mix_dist = self._normal_dist(mix, **assay_params)
        pool_dist = self._normal_dist(pool, **assay_params)

        # Step 2: Select random sets of parameters from the curve_fit model
        if use_method == 'curve':
            curve_fit_param_indices = np.random.randint(0, self.out.shape[0]-1, size=n_iters)
            curve_fit_params = self.out[curve_fit_param_indices, :].T
            # Step 3: Transform the Straight/Mix/Pool values (seconds) to factor level, 
            # based on the curve_fit_params selected in step 2
            # _fl == "_factor_level"
            straight_fl = self._fwd_func(straight_dist, *curve_fit_params)
            mix_fl = self._fwd_func(mix_dist, *curve_fit_params)
            pool_fl = self._fwd_func(pool_dist, *curve_fit_params)
        elif use_method == 'linear':
            straight_fl = self._fwd_func_linear(straight_dist)
            mix_fl = self._fwd_func_linear(mix_dist)
            pool_fl = self._fwd_func_linear(pool_dist)
        else:
            raise ModelException("unknown method: %s" % use_method)

        # Step 4: Calculate the target factor level (Target_fl)
        # And then the perecent_corrected
        target_fl = (straight_fl + pool_fl) / 2.0
        percent_corrected_dist = self._percent_corrected(straight_fl, mix_fl, pool_fl, target_fl)

        return percent_corrected_dist

    def _normal_dist(self, value, sigma, n_iters):
        return np.random.normal(loc=value, scale=sigma, size=n_iters)