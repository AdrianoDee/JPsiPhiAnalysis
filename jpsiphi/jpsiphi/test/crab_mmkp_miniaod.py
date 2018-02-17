import sys
import os

jsonFile="Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt"

from WMCore.Configuration import Configuration
config = Configuration()

#print("Test = " + str(skipevt))

datasetbase = '/Charmonium' # '/Muonia' #

sites = ['T2_AT_Vienna', 'T2_BE_IIHE', 'T2_BE_UCL', 'T2_BR_SPRACE', 'T2_BR_UERJ',
 'T2_CH_CERN', 'T2_CH_CERN_AI', 'T2_CH_CERN_HLT',
 'T2_CH_CSCS', 'T2_CH_CSCS_HPC', 'T2_CN_Beijing', 'T2_DE_DESY', 'T2_DE_RWTH',
 'T2_EE_Estonia', 'T2_ES_CIEMAT', 'T2_ES_IFCA', 'T2_FI_HIP', 'T2_FR_CCIN2P3',
 'T2_FR_GRIF_IRFU', 'T2_FR_GRIF_LLR', 'T2_FR_IPHC', 'T2_GR_Ioannina', 'T2_HU_Budapest',
 'T2_IN_TIFR', 'T2_IT_Bari', 'T2_IT_Legnaro', 'T2_IT_Pisa', 'T2_IT_Rome', 'T2_KR_KISTI',
 'T2_KR_KNU', 'T2_MY_UPM_BIRUNI', 'T2_PK_NCP', 'T2_PL_Swierk', 'T2_PL_Warsaw',
 'T2_PT_NCG_Lisbon', 'T2_RU_IHEP', 'T2_RU_INR', 'T2_RU_ITEP', 'T2_RU_JINR', 'T2_RU_PNPI',
 'T2_RU_SINP', 'T2_TH_CUNSTDA', 'T2_TR_METU', 'T2_TW_NCHC', 'T2_UA_KIPT', 'T2_UK_London_Brunel',
 'T2_UK_London_IC', 'T2_UK_SGrid_Bristol','T2_UK_SGrid_RALPP', 'T2_US_Caltech', 'T2_US_Florida',
 'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD', 'T2_US_Vanderbilt',
 'T2_US_Wisconsin',

 'T3_CH_PSI', 'T3_CN_PKU', 'T3_CO_Uniandes',
 'T3_ES_Oviedo', 'T3_FR_IPNL', 'T3_GR_IASA', 'T3_HU_Debrecen',
 'T3_IN_PUHEP', 'T3_IN_TIFRCloud','T3_IT_*', 'T3_KR_*', 'T3_MX_*', 'T3_RU_*',
 #'T3_TW*', #Taiwan gives some troubles with file fetching
 'T3_UK_*', 'T3_US_Baylor','T3_US_Colorado', 'T3_US_Cornell',
 'T3_US_FIT', 'T3_US_FIU', 'T3_US_FNALLPC', 'T3_US_FSU', 'T3_US_J*',
 'T3_US_Kansas', 'T3_US_MIT', 'T3_US_N*', 'T3_US_O*', 'T3_US_P*',
 'T3_US_R*', 'T3_US_S*', 'T3_US_T*', 'T3_US_UCD', 'T3_US_UCR', 'T3_US_UCSB',
 'T3_US_UMD']

datasetnames = {

"F" :  datasetbase + '/Run2017F-PromptReco-v1/MINIAOD',
"B1" : datasetbase + '/Run2017B-PromptReco-v1/MINIAOD',
"B2" : datasetbase + '/Run2017B-PromptReco-v2/MINIAOD',
"C1" : datasetbase + '/Run2017C-PromptReco-v1/MINIAOD',
"C2" : datasetbase + '/Run2017C-PromptReco-v2/MINIAOD',
"C3" : datasetbase + '/Run2017C-PromptReco-v3/MINIAOD',
"D" : datasetbase + '/Run2017D-PromptReco-v1/MINIAOD',
"E" : datasetbase + '/Run2017E-PromptReco-v1/MINIAOD'
}


runNumber = [
#'274094-274240',
'',
]

run = 'F'

datasetName = datasetnames[run]
runNum = runNumber[0]
#lumi = jsonfile[jNum]
lumi = jsonFile
#HLT = HLTPath[0]

import datetime
timestamp = datetime.datetime.now().strftime("_%Y%m%d_%H%M%S")

dataset = filter(None, datasetName.split('/'))

jobdir = 'miniaod_X4140_MMKPI_' + run

if not os.path.exists(jobdir):
    os.makedirs(jobdir)

config.section_('General')
config.General.transferOutputs  = True
config.General.workArea         = jobdir
#config.General.requestName     = 'JetHT_Run2015D_PromptReco_v4_RECO'+timestamp
#config.General.requestName             = dataset[0]+'_'+dataset[1]+'_'+dataset[2]+'_'+runNum+'_'+HLT+timestamp
config.General.requestName      = 'mmpik_miniaod_' + dataset[0]+'_'+dataset[1]+'_'+dataset[2]+'_'+runNum+'_'+timestamp #+'_split_'+ jsonFile.split('_')[-1].split('.')[0]
config.General.transferLogs     = False

config.section_('JobType')
config.JobType.psetName         = '/lustre/home/adrianodif/CMSSW_9_2_13/src/mmkk/mmkk/test/run-mumupik-miniaod.py'
config.JobType.pluginName       = 'Analysis'
config.JobType.maxMemoryMB      = 2500
config.JobType.maxJobRuntimeMin = 2750
config.JobType.allowUndistributedCMSSW = True
#config.JobType.inputFiles      = ['Run_time.txt']
#config.JobType.outputFiles     = ['EventList.txt','EventListClean.txt']

config.section_('Data')
config.Data.inputDataset        = datasetName
config.Data.inputDBS            = 'global'
config.Data.totalUnits          = -1
config.Data.unitsPerJob         = 6
config.Data.splitting           = 'LumiBased'
config.Data.runRange            = runNum
config.Data.lumiMask            = lumi
config.Data.outLFNDirBase       = '/store/user/adiflori/'
config.Data.publication         = False
config.Data.ignoreLocality      = True


config.section_('Site')
#config.Site.storageSite        = 'T2_CH_CERN'
config.Site.storageSite         = 'T2_IT_Bari'
#config.Site.blacklist          = ['T2_IN_TIFR','T2_US_Vanderbilt']
config.Site.whitelist           = sites
config.Site.blacklist           = ['T1*', 'T3_US_UMiss']
