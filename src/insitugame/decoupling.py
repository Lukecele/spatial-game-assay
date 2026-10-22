"""
decoupling.py
=============
Decouples macroenvironmental physical gradients (hypoxia, vasculature)
from local spot proliferation to obtain purified evolutionary game fitness.
"""

import numpy as np
from sklearn.linear_model import LinearRegression


class EnvironmentalDecoupler:
    def __init__(self):
        self.model = LinearRegression()
        self.r2_env = 0.0
        
    def fit_transform(self, w_raw: np.ndarray, z_hyp: np.ndarray, z_vasc: np.ndarray) -> np.ndarray:
        X_env = np.column_stack([z_hyp, z_vasc])
        self.model.fit(X_env, w_raw)
        self.r2_env = float(self.model.score(X_env, w_raw))
        w_purified = w_raw - self.model.predict(X_env)
        return w_purified
