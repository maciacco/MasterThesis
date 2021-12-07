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

def expected_signal(cent_class, ct_range, eff, n_events, cent_counts, cent_edges):
    lambda_yield_list = [2.35e-4, 2.03e-4, 6.58e-5] # TODO: update centrality-differential yields
    correction = expo(ct_range[0])- expo(ct_range[1]) #selecting the correct ct bin
    correction *= eff
    correction *= 0.639 # branching ratio
    correction *= 2 # antimatter + matter
    cent_end_bin = [5., 10., 50.]
    for cent_bin, lambda_yield in zip(cent_end_bin, lambda_yield_list):
        if cent_bin==cent_class[1]:
            return lambda_yield*correction*n_events

    # 0-90% production in Pb-Pb at 5.02 TeV
    lambda_yield_list_0_90 = [32.30,26.48,20.60,14.36,9.69,6.07,3.57,1.84,0.89,0.33]
    cent_list = [[0,5],[5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]]
    
    # compute yield
    lambda_yield_0_90 = 0.
    cent_bin_centers = (cent_edges[:-1]+cent_edges[1:])/2
    for i_cent_bins in range(len(cent_list)):
        # print(f"index = {i_cent_bins}")
        cent_range_map = np.logical_and(cent_bin_centers > cent_list[i_cent_bins][0], cent_bin_centers < cent_list[i_cent_bins][1])
        counts_cent_range = cent_counts[cent_range_map]
        counts_cent_range_sum = np.sum(counts_cent_range)
        lambda_yield_0_90 += counts_cent_range_sum*lambda_yield_list_0_90[i_cent_bins]
        
    # expected signal for 0-90% centrality
    return correction*lambda_yield_0_90