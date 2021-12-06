#!/usr/bin/env python3
import os
import pickle
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import uproot
import yaml


def ndarray2roo(ndarray, var):
    if isinstance(ndarray, ROOT.RooDataSet):
        print('Already a RooDataSet')
        return ndarray

    assert isinstance(ndarray, np.ndarray), 'Did not receive NumPy array'
    assert len(ndarray.shape) == 1, 'Can only handle 1d array'

    name = var.GetName()
    x = np.zeros(1, dtype=np.float64)

    tree = ROOT.TTree('tree', 'tree')
    tree.Branch(f'{name}', x, f'{name}/D')

    for i in ndarray:
        x[0] = i
        tree.Fill()

    array_roo = ROOT.RooDataSet(
        'data', 'dataset from tree', tree, ROOT.RooArgSet(var))
    return array_roo


def significance_error(signal, background):
    num_s_prop = 0.5*signal+background
    num_b_prop = 0.5*signal
    den_prop = np.sqrt((signal+background)*(signal+background)*(signal+background))
    return 1/den_prop*np.sqrt(signal*num_s_prop*num_s_prop+background*num_b_prop*num_b_prop)

def expo(x):
    return np.exp(-x / (252 * 0.029979245800)) #hyp tau taken from lifetime analysis

def expected_signal(cent_class, ct_range, eff, n_events):
    lambda_yield_list = [2.35e-4, 2.03e-4, 6.58e-5] # TODO: update centrality-differential yields
    correction = expo(ct_range[0])- expo(ct_range[1]) #selecting the correct ct bin
    correction *= eff
    cent_end_bin = [5., 10., 50.]
    for cent_bin, lambda_yield in zip(cent_end_bin, lambda_yield_list):
        if cent_bin==cent_class[1]:
            return lambda_yield*correction*n_events

    # expected signal for 0-90% centrality
    lambda_yield_0_90 = 1.55e8/55
    return correction*lambda_yield_0_90