#perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/DYJets_mumuFilter_test.txt isMC=True --prodSpace /data/cmszfs1/user/revering/ZmmSim --batch 10 --jobname ZmmSim_testBarrelCheck
perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/DYJets_mumuFilter_combined.txt isMC=True --prodSpace /local/cms/user/auerb029/ZmmSim --batch 10 --jobname ZmmSim_testBarrelCheck
perl condor_filelist.perl python/ConfFile_MuAnalyzer_cfg.py datafiles/DYJets_mumuFilter_part2.txt isMC=True --prodSpace /local/cms/user/auerb029/ZmmSim --batch 10 --jobname ZmmSim_testBarrelCheck_2
