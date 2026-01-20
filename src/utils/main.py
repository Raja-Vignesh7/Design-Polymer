import pandas as pd
import numpy as np
import warnings
from rdkit import Chem
from rdkit.Chem import AllChem
import warnings
import sys
import os
from rdkit import rdBase
from mordred import Calculator, descriptors
from catboost import  CatBoostRegressor
import lightgbm as lgb
from lightgbm import  LGBMRegressor
from sklearn.ensemble import VotingRegressor, AdaBoostRegressor
from xgboost import XGBRegressor
from catboost import CatBoostRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, explained_variance_score
from sklearn.model_selection import RandomizedSearchCV
import matplotlib.pyplot as plt
import numpy as np
from SMILE_handler import get_smile_features
import warnings
import joblib
warnings.filterwarnings('ignore')


class Models:
    def __init__(self):
        self.tg_model = joblib.load('src\models\Tg_model.joblib')
        self.Tc_model = joblib.load('src\models\Tc_model.joblib')
        self.Rg_model = joblib.load('src\models\Rg_model.joblib')
        self.FFV_model = joblib.load('src\models\FFV_model.joblib')
        self.Density_model = joblib.load('src\models\Density_model.joblib')
        
    def predict_properties(self, smile):
        features = get_smile_features(smile)
        tg_pred = self.tg_model.predict([features])[0]
        Tc_pred = self.Tc_model.predict([features])[0]
        Rg_pred = self.Rg_model.predict([features])[0]
        FFV_pred = self.FFV_model.predict([features])[0]
        Density_pred = self.Density_model.predict([features])[0]
        
        return {
            "Tg": tg_pred,
            "Tc": Tc_pred,
            "Rg": Rg_pred,
            "FFV": FFV_pred,
            "Density": Density_pred
        }
        
# ip = "*CC(*)c1ccccc1C(=O)OCCCCCC"  # Benzene
# model = Models()
# print(model.predict_properties(ip))