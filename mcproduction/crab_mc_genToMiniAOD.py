from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import datetime, time
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d-%H-%M')

step = 'GEN-MINIAODSIM'
nEvents = 400000
job_label = 'generic'
myrun= 'step0-GS-generic_cfg.py'
myname=step+job_label

config.General.requestName = step+'-'+job_label+'-'+st
config.General.transferOutputs = True
config.General.transferLogs = False
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = myrun
config.JobType.inputFiles = ['step0-GS-generic_cfg.py','step1-DR-generic_cfg.py','step2-DR-generic_cfg.py','step3-MiniAOD-generic_cfg.py','Ponia.tar.gz']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.eventsPerLumi=10000
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 3000
config.JobType.scriptExe = 'GEN-MiniAOD-generic.sh'#'GEN-MiniAOD-Xb2chib1pipi_10p5.sh'
#config.JobType.scriptArgs = ['=Xb2chib1pipi','=10p5','=100000']
config.JobType.outputFiles = ['step3-MiniAOD-generic.root','Rootuple-GEN-MiniAOD-2017-generic.root','Rootuple-MiniAOD-2017-generic.root']

config.Data.outputPrimaryDataset = myname
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = nEvents # the number of events here must match the number of events in the externalLHEProducer
NJOBS = 2000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.outputDatasetTag = 'RunIISummer17PrePremix-MC_v2_94X_mc2017_realistic_v11'+step
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.publication = False

config.Site.storageSite = 'T2_IT_Bari'
#config.Site.whitelist = ['T2_IT_*','T2_CH_CERN','T2_IT_Bari','T2_IT_Legnaro','T2_IT_Pisa','T2_IT_Rome']
