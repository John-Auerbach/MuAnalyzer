To initialize:

cmsrel CMSSW_10_6_20

cd CMSSW_10_6_20/src/

mkdir DarkPhoton

cd DarkPhoton/

git clone (this repository)


To compile:

cmsenv
scram b

To run:
submitScripts/submitDYMC.sh
