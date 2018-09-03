from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import datetime, time
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M%S')

sites = ['T2_AT_Vienna', 'T2_BE_IIHE', 'T2_BE_UCL', 'T2_BR_SPRACE', 'T2_BR_UERJ',

 'T2_DE_DESY', 'T2_DE_RWTH',
 'T2_EE_Estonia', 'T2_ES_CIEMAT', 'T2_ES_IFCA', 'T2_FI_HIP', 'T2_FR_CCIN2P3',
 'T2_FR_GRIF_IRFU', 'T2_FR_GRIF_LLR', 'T2_FR_IPHC', 'T2_GR_Ioannina', 'T2_HU_Budapest',
 'T2_IN_TIFR', 'T2_IT_Bari', 'T2_IT_Legnaro', 'T2_IT_Pisa', 'T2_IT_Rome', 'T2_KR_KISTI',
  'T2_PK_NCP', 'T2_PL_Swierk', 'T2_PL_Warsaw',
 'T2_PT_NCG_Lisbon', 'T2_RU_IHEP', 'T2_RU_INR', 'T2_RU_ITEP',
 'T2_RU_SINP','T2_TR_METU', 'T2_TW_NCHC', 'T2_UA_KIPT', 'T2_UK_London_Brunel',
 'T2_UK_London_IC', 'T2_UK_SGrid_Bristol','T2_UK_SGrid_RALPP', #'T2_US_Caltech', 'T2_US_Florida',
 'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD', 'T2_US_Vanderbilt',
 'T2_US_Wisconsin',

# 'T3_CH_PSI', 'T3_CN_PKU', 'T3_CO_Uniandes',
 'T3_ES_Oviedo', 'T3_GR_IASA', 'T3_HU_Debrecen',
# 'T3_IN_PUHEP', 'T3_IN_TIFRCloud',
 'T3_IT_*',
 #'T3_KR_*', 'T3_MX_*', 'T3_RU_*',
 #'T3_TW*', #Taiwan gives some troubles with file fetching
 'T3_UK_*'#, 'T3_US_Baylor','T3_US_Colorado', 'T3_US_Cornell',
# 'T3_US_FIT', 'T3_US_FIU', 'T3_US_FNALLPC', 'T3_US_FSU', 'T3_US_J*',
# 'T3_US_Kansas', 'T3_US_MIT', 'T3_US_N*', 'T3_US_O*', 'T3_US_P*',
# 'T3_US_R*', 'T3_US_S*', 'T3_US_T*', 'T3_US_UCD', 'T3_US_UCR', 'T3_US_UCSB',
# 'T3_US_UMD']
]

NJOBS = 2500
EVTPERJOB = 200000
step = 'GEN-MINIAODSIM'
nEvents = NJOBS*EVTPERJOB
ptHat = 50
genname = 'BBbar_JpsiFilter_SoftQCD'
job_label = 'BBbar_JpsiFilter_SoftQCD'
myrun= job_label + '.py'
myname=step+job_label

step1File = genname + '_GEN_SIM.root'
step2File = genname + '_DIGI_HLT_RAW_L1_PU40.root'
step3Mini = genname + '_MINIAODSIM_PU40.root'
step3File = genname + '_RECOSIM_PU40.root'
config.General.requestName = step+'_'+job_label+'_'+st
config.General.transferOutputs = True
config.General.transferLogs = False
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = myrun
config.JobType.inputFiles = [myrun,'step2_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_PU.py','step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_PU.py']#,'Ponia.tar.gz']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.eventsPerLumi=10000
#config.JobType.numCores = 1
config.JobType.maxMemoryMB = 3000
config.JobType.scriptExe = 'mcproduction.sh'#'GEN-MiniAOD-Xb2chib1pipi_10p5.sh'
#config.JobType.scriptArgs = ['=Xb2chib1pipi','=10p5','=100000']
config.JobType.outputFiles = [step1File,step3Mini]

config.Data.outputPrimaryDataset = myname
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = EVTPERJOB#nEvents # the number of events here must match the number of events in the externalLHEProducer

config.Data.totalUnits = nEvents
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.outputDatasetTag = 'RunIISummer17PrePremix-MC_v2_94X_mc2017_realistic_v11'+step
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.publication = True
config.Data.inputDBS = 'phys03'
config.Data.publishDBS = 'phys03'

config.Site.storageSite = 'T2_IT_Bari'
config.Site.blacklist = ["T2_BR_*","T2_IN_*","T3_IN_*"]
#config.Site.whitelist = sites #['T2_IT_*','T2_CH_*','T2_GE_*','T2_FR_*','T2_ES_*','T2_UK_*']
