import os
import pandas as pd 
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/table_data.tsv")
TABLES = pd.read_csv(DATA_PATH, sep="\t", index_col="table")

class MixingConverter(object):
    """docstring for MixingConverter"""
    def __init__(self, table_name="vitk", method="curve_fit"):
        super(MixingConverter, self).__init__()
        self.table_name = table_name
        self.method = method

        self.table = TABLES.loc[table_name]
        self.PT_POOL_SECONDS = 13.6
        self.PTT_POOL_SECONDS = 28.6
        
        if self.method == "curve_fit":
            def fwd_func(x, a, c, d):
                return a*np.exp(-c*x)+d

            def rev_func(y, a, c, d):
                return np.log((y-d)/a)/-c

            self.ptt_params, _ = curve_fit(fwd_func, self.table.PTT, self.table.factor_level, p0=(40000, 0.1, 20))
            self.pt_params, _ = curve_fit(fwd_func, self.table.PT, self.table.factor_level, p0=(40000, 0.1, 20))

            self.pt_to_level   = lambda seconds: fwd_func(seconds, *self.pt_params)
            self.ptt_to_level  = lambda seconds: fwd_func(seconds, *self.ptt_params)
            self.level_to_pt   = lambda level: rev_func(level, *self.pt_params)
            self.level_to_ptt  = lambda level: rev_func(level, *self.ptt_params)

        elif self.method == "linear":
            ptt_spl = UnivariateSpline(self.table.PTT, self.table.factor_level, k=1,s=0)
            pt_spl = UnivariateSpline(self.table.PT, self.table.factor_level, k=1,s=0)
            # It is not totally clear why these values need to be reversed for the spline to be fit
            ptt_spl_rev = UnivariateSpline(self.table.factor_level[::-1], self.table.PTT[::-1], k=1,s=0)
            pt_spl_rev = UnivariateSpline(self.table.factor_level[::-1], self.table.PT[::-1], k=1,s=0)

            self.pt_to_level   = lambda seconds: pt_spl(seconds)
            self.ptt_to_level  = lambda seconds: ptt_spl(seconds)
            self.level_to_pt   = lambda level: pt_spl_rev(level)
            self.level_to_ptt  = lambda level: ptt_spl_rev(level)

        else:
            raise Exception("Invalid `method` for class MixingConverter")

    def pt_percent_corrected(self, pt_straight, pt_mix, pt_pool=None):
        if np.isnan(pt_mix) or np.isnan(pt_straight):
            return None
        pt_pool = pt_pool or self.PT_POOL_SECONDS
        levels = self.pt_to_level(np.array([pt_straight, pt_mix, pt_pool]))
        return self._percent_corrected(*levels)

    def ptt_percent_corrected(self, ptt_straight, ptt_mix, ptt_pool=None):
        if np.isnan(ptt_mix) or np.isnan(ptt_straight):
            return None
        ptt_pool = ptt_pool or self.PTT_POOL_SECONDS
        levels = self.ptt_to_level(np.array([ptt_straight, ptt_mix, ptt_pool]))
        return self._percent_corrected(*levels)

    def oldway_combined_percent_corrected(self, pt_straight, ptt_straight, pt_mixed, ptt_mixed):

        # Step 1: Convert seconds to level
        pt_straight_level, pt_mixed_level = self.pt_to_level([pt_straight, pt_mixed])
        ptt_straight_level, ptt_mixed_level = self.ptt_to_level([ptt_straight, ptt_mixed])

        # Step 2: Take minimum of both PT and PTT factor level
        # (per LMR call manual)
        factor_level = min(pt_straight_level, ptt_straight_level)

        # Step 3: Determine target level, assuming 100% pool level
        # (per LMR call manual)
        target_level = self._target_level(factor_level, 100)
        
        # Step 4: Calculate percentage correction for each, if relevant
        if pt_straight_level < target_level:
            pt_percent_correct = self._percent_corrected(pt_straight_level, pt_mixed_level, None, target_level)
        else:
            pt_percent_correct = None

        if ptt_straight_level < target_level:
            ptt_percent_correct = self._percent_corrected(ptt_straight_level, ptt_mixed_level, None, target_level)
        else:
            ptt_percent_correct = None

        return pt_percent_correct, ptt_percent_correct
    
    def _percent_corrected(self, straight, mix, pool=100, target=None):
        target = target or self._target_level(straight, pool)
        try:
            return 100. * (mix - straight) / (target - straight)
        except ZeroDivisionError:
            return 100
    
    def _target_level(self, initial, pool):
        return (initial + pool) / 2.
    
    def _plot_table(self, column, ax):
        import matplotlib.pyplot as plt

        if column == "PT":
            func = self.pt_to_level
        elif column == "PTT":
            func = self.ptt_to_level
        else:
            raise Exception("Invalid column type")
        ax.scatter(self.table[column], self.table.factor_level, label=self.table_name)
        xp = np.linspace(0, 250, 250)
        ax.plot(xp, func(xp), '-', lw=1)
        ax.set_ylim([0,170])
        ax.set_xlim([0,100])
        ax.set_title(column)
        ax.set_xlabel("Seconds")
        ax.set_ylabel("% Factor Present")
        
    def plot_ptt_table(self, ax=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(6,5))
        self._plot_table("PTT", ax)

    def plot_pt_table(self, ax=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(6,5))
        self._plot_table("PT", ax)
    
    def plot_tables(self):
        fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(12,5))
        self.plot_pt_table(ax=axes[0])
        self.plot_ptt_table(ax=axes[1])
    
    def pp_table(self, row):
        out = pd.DataFrame(index=["Patient", "Mix", "Pool"], 
                           columns=["PT", "PTT", "PT_level", "PTT_level"])
        out.loc[:,"PT"] = [row["PROPT"], row["PROMIX"], row["PRONOR"]]
        out.loc[:,"PTT"] = [row["PTTPT"], row["PTTMIX"], row["PTTNOR"]]
        try:
            out.loc[:,"PT_level"] = np.round(self.pt_to_level(out.PT),0).astype(int)
        except:
            out.loc[:,"PT_level"] = None
        try:
            out.loc[:,"PTT_level"] = np.round(self.ptt_to_level(out.PTT),0).astype(int)
        except:
            out.loc[:,"PTT_level"] = None

        return out