"""Regionalised parameterisation theta = gA(A, c) (Feng et al. 2022, Sec. 2.1).

Copied from delta_gwb/dpl.py. One change in how it is called, not in the
network: the dynamic input is the cell's day-of-year climate normals, a fixed
sequence, so the "static" parameters really are the same for every period
(finding L5). Before, they were re-inferred from whichever window was being
simulated, including the test window.
"""
import torch
import torch.nn as nn


class ParamNet(nn.Module):
    def __init__(self, n_attr, n_param, hidden=256, layers=2, dropout=0.0, n_dyn=3):
        super().__init__()
        self.lstm = nn.LSTM(input_size=n_dyn + n_attr, hidden_size=hidden, num_layers=layers,
                            dropout=dropout if layers > 1 else 0.0, batch_first=True)
        self.head = nn.Linear(hidden, n_param)

    def forward(self, dyn, A):
        """dyn (B, T, n_dyn), A (B, n_attr) -> raw parameters (B, n_param)."""
        sta = A.unsqueeze(1).expand(-1, dyn.shape[1], -1)
        _, (h, _) = self.lstm(torch.cat([dyn, sta], dim=-1))
        return self.head(h[-1])
