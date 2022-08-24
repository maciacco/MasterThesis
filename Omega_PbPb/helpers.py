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

# from DmesonAnalysers/DmesonAnalysis
def GetPromptFDYieldsAnalyticMinimisation(effPromptList, effFDList, rawYieldList, effPromptUncList, effFDUncList,
                                          rawYieldUncList, corr=True, precision=1.e-8, nMaxIter=100):
    '''
    Method to retrieve prompt and FD corrected yields with an analytic system minimisation
    Parameters
    ----------
    - effPromptList: list of efficiencies for prompt D
    - effFDList: list of efficiencies for FD D
    - rawYieldList: list of raw yields
    - effPromptUncList: list of uncertainties on efficiencies for prompt D
    - effFDUncList: list of uncertainties on efficiencies for FD D
    - rawYieldUncList: list of uncertainties on raw yields
    - corr (bool, optional): whether to compute the correlation
    - precision (float, optional): target precision for minimisation procedure
    - nMaxIter (int, optional): max number of iterations for minimisation procedure
    Returns
    ----------
    - mCorrYield (numpy matrix): corrected yields (Nprompt, NFD)
    - mCovariance (numpy matrix): covariance matrix for corrected yields
    - redChiSquare (float): reduced chi square
    - dicOfMatrices (dictionary): dictionary with all matrices used in minimisation procedure
    '''
    nCutSets = len(effPromptList)

    mRawYield = np.zeros(shape=(nCutSets, 1))
    mEff = np.zeros(shape=(nCutSets, 2))
    mCovSets = np.zeros(shape=(nCutSets, nCutSets))
    mCorrSets = np.zeros(shape=(nCutSets, nCutSets))
    mWeights = np.zeros(shape=(nCutSets, nCutSets))

    mCorrYield = np.zeros(shape=(2, 1))
    mCorrYieldOld = np.zeros(shape=(2, 1))
    mCovariance = np.zeros(shape=(2, 2))
    mRes = np.zeros(shape=(nCutSets, 1))

    for iCutSet, (rawYield, effPrompt, effFD) in enumerate(zip(rawYieldList, effPromptList, effFDList)):
        mRawYield.itemset(iCutSet, rawYield)
        mEff.itemset((iCutSet, 0), effPrompt)
        mEff.itemset((iCutSet, 1), effFD)

    mRawYield = np.matrix(mRawYield)
    mEff = np.matrix(mEff)

    for iIter in range(nMaxIter):
        if iIter == 0:
            mCorrYield.itemset(0, 0)
            mCorrYield.itemset(1, 0)
        for iCutSetRow, (rawYieldUncRow, effPromptUncRow, effFDUncRow) in enumerate(\
            zip(rawYieldUncList, effPromptUncList, effFDUncList)):
            for iCutSetCol, (rawYieldUncCol, effPromptUncCol, effFDUncCol) in enumerate(\
                zip(rawYieldUncList, effPromptUncList, effFDUncList)):
                uncRow = np.sqrt(rawYieldUncRow**2 + effPromptUncRow**2 *
                                 mCorrYield.item(0)**2 + effFDUncRow**2 * mCorrYield.item(1)**2)
                uncCol = np.sqrt(rawYieldUncCol**2 + effPromptUncCol**2 *
                                 mCorrYield.item(0)**2 + effFDUncCol**2 * mCorrYield.item(1)**2)
                if corr and uncRow > 0 and uncCol > 0:
                    if uncRow < uncCol:
                        rho = uncRow / uncCol
                    else:
                        rho = uncCol / uncRow
                else:
                    if iCutSetRow == iCutSetCol:
                        rho = 1.
                    else:
                        rho = 0.
                covRowCol = rho * uncRow * uncCol
                mCovSets.itemset((iCutSetRow, iCutSetCol), covRowCol)
                mCorrSets.itemset((iCutSetRow, iCutSetCol), rho)

        mCovSets = np.matrix(mCovSets)
        mWeights = np.linalg.inv(np.linalg.cholesky(mCovSets))
        mWeights = mWeights.T * mWeights
        mEffT = mEff.T

        mCovariance = (mEffT * mWeights) * mEff
        mCovariance = np.linalg.inv(np.linalg.cholesky(mCovariance))
        mCovariance = mCovariance.T * mCovariance

        mCorrYield = mCovariance * (mEffT * mWeights) * mRawYield
        mRes = mEff * mCorrYield - mRawYield
        mResT = np.transpose(mRes)

        if (mCorrYield.item(0)-mCorrYieldOld.item(0)) / mCorrYield.item(0) < precision and \
            (mCorrYield.item(1)-mCorrYieldOld.item(1)) / mCorrYield.item(1) < precision:
            break

        mCorrYieldOld = np.copy(mCorrYield)

    #reduced chi2
    redChiSquare = mResT * mWeights * mRes / (nCutSets - 2)
    #dictionary with matrices used in minimisation procedure
    dicOfMatrices = {'covMatrix':mCovSets, 'weightMatrix':mWeights, 'corrMatrix':mCorrSets}

    return mCorrYield, mCovariance, float(redChiSquare), dicOfMatrices

