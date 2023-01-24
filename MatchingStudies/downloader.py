import os

data_type = ['mc']
data_path = '/data/mciacco/MatchingStudies'
train_id_mc = '1058_20230109-1304'
train_id_dat = ['1434_20230109-1105', '1433_20230109-1105']
mother = 'K0s'
period = [['q', 1], ['r', 2]]

if not os.path.exists(f'{data_path}'):
    os.makedirs(f'{data_path}')

for dat in data_type:
    if not os.path.exists(f'{data_path}/{dat}'):
        os.makedirs(f'{data_path}/{dat}')
    if dat == 'mc':
        for p in period:
            for child in range(1,4):
                if child == 2 and p[0] == 'q':
                    continue
                os.system(f'alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGLF/LF_PbPb_MC_AOD/{train_id_mc}_child_{child}/merge_runlist_{p[1]}/AnalysisResults.root file:{data_path}/{dat}')
                os.system(f'mv {data_path}/{dat}/AnalysisResults.root {data_path}/{dat}/AnalysisResults_{child}.root')
            os.system(f'alihadd {data_path}/{dat}/AnalysisResults.root {data_path}/{dat}/AnalysisResults_1.root {data_path}/{dat}/AnalysisResults_2.root {data_path}/{dat}/AnalysisResults_3.root')
            os.system(f'rm {data_path}/{dat}/AnalysisResults_*')
            os.system(f'mv {data_path}/{dat}/AnalysisResults.root {data_path}/{dat}/AnalysisResults_{p[0]}_{mother}.root')
    else:
        for t, p in zip(train_id_dat, period):
            os.system(f'alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGLF/LF_PbPb_AOD/{t}/merge/AnalysisResults.root file:{data_path}/{dat}')
            os.system(f'mv {data_path}/{dat}/AnalysisResults.root {data_path}/{dat}/AnalysisResults_{p[0]}_{mother}.root')

