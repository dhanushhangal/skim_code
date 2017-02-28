from WMCore.Configuration import Configuration
config = Configuration()

#subScript = "unzipVectorTree.C"
subScript = "mini_skim.C"
#subScript = "HT_Analyzer_All_JFFCorr2.C"

config.section_("General")
config.General.requestName = 'PbPb_Pythia6Hydjet_unzippedSkim_forJetPtClosures_testClosure_24Feb2017'
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'runScript_unzip.sh'
#config.JobType.scriptExe = 'runScript.sh'
config.JobType.scriptArgs = ['script='+subScript]
config.JobType.inputFiles = ['FrameworkJobReport.xml',subScript]
config.JobType.outputFiles = ['unzippedSkim_new.root']
config.JobType.maxJobRuntimeMin = 600

config.section_("Data")
#config.Data.userInputFiles = open('PbPb_MC_Pythia6Hydjet_5TeV.txt').readlines()
#config.Data.userInputFiles = open('PbPb_Hydjet_MC_Jan4Skim_newjetCorr.txt').readlines() 
config.Data.userInputFiles = open('Pythia6+Hydjet.txt').readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'PbPb_MC_Histograms'
#config.Data.outLFNDirBase = '/afs/cern.ch/work/d/dhangal/ncs_skims/'+config.General.requestName
config.Data.outLFNDirBase = '/store/group/phys_heavyions/dhangal/ncs_skims'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = 'T2_CH_CERN'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

