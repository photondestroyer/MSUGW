"""Residual-ET (groundwater-buffering) pipeline shared by the XGBoost and the
differentiable LSTM-HBV models. See RESIDUAL_ET_MODELS.md at the repo root.

Stages:  features (prepare)  ->  train_xgb / train_dlstm  ->  buffering
Sources under merged_datasets/ are opened read-only; everything written goes
to the store directory (config.STORE).
"""
