import sys
import os
import math

# MEM Base
MEM_BASE = os.environ['CMSSW_BASE']+ "/src/ttH_Htautau_MEM_Analysis/MEM"

# Math constants
PI          = math.pi

# LHC beam energy (Gev)
sqrtS = 13000



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Input files
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

InputType = "PyRun2"

InputFileList = ["ttH_Htautau_MEM_Analysis/MEM/EventDump_ttH.py"]



#treeName = "HTauTauTree_2lSS"


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Events
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
maxNbrOfEventsToRead = 2
skipEvents = 0


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Random Number Generator
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
flagSameRNG = True

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Computation / Run
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

flagTFLepTau = True 
flagTFHadTau = True
flagTFFake = True
flagTFMET    = True
flagTFJet1   = True
flagTFJet2   = True
flagTFBJet_leptop = True
flagTFBJet_hadtop = True
flagTF_fakelep = True
flagTF_fakeleptau = True
flagTFTop    = True
flagJac      = True
flagWME      = True
#  Verbose Modes
NoVerbose        = 0
ResultLevel      = 1 
IntegrationLevel = 2 
IntegrandLevel   = 3 
# Seclect verbose mode
verbose      = NoVerbose
#logFileName = MEM_BASE + "/python/small.log"



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MC Integration
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
nbrOfDimttH    = 5
nbrOfDimttZ    = 5
nbrOfDimttbar_DL = 5
nbrOfDimttZ_Zll    = 3


nbrOfDimttH_miss    = 7
nbrOfDimttZ_miss    = 7
nbrOfDimttZ_Zll_miss    = 5


runttZ_integration = True
runttZ_Zll_integration = True
runttbar_DL_fakelep_integration = True

run_missing_jet_integration = True
force_missing_jet_integration = False
force_missing_jet_integration_ifnoperm = True



nbrOfPoints_ttH = 10000
#nbrOfPoints_ttH = 15
nbrOfPoints_ttZ = 10000
#nbrOfPoints_ttZ = 10
nbrOfPoints_ttbar_DL = 10000
#nbrOfPoints_ttbar_DL = 10
nbrOfPoints_ttZ_Zll = 5000
#nbrOfPoints_ttZ_Zll = 5


nbrOfPoints_ttH_miss = 20000
#nbrOfPoints_ttH_miss = 20
nbrOfPoints_ttZ_miss = 20000
#nbrOfPoints_ttZ_miss = 20
nbrOfPoints_ttZ_Zll_miss = 10000
#nbrOfPoints_ttZ_Zll_miss = 10

nbrOfPermut_per_jet = 4

# Integration boundaries 
phiNu_tlep  = [ 0, 2*PI ]
cosThetaNu_tlep = [-1.,1.]
phiNu_ttau  = [ 0, 2*PI ]
cosThetaNu_ttau = [-1.,1.]
 
use_pT_TFJet = True
CI_TFJet = 0.95 #Use 95% CI from the jet TF
use_top_compatibility_check = True

include_hadrecoil = True

#Missing jet parameters
eta_acceptance = 2.4
jet_radius = 0.4
dR_veto_jet_lep = 0.4
rel_iso_lep = 0.4
pT_cut = 25.

phi_missing_jet = [ 0, 2*PI ]
cosTheta_missing_jet = [-1.,1.]



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MadGraph 5
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
MGParamCardFile =  MEM_BASE+"/bin/param_card.dat"


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# LHAPDF
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
LHAPDFFile = "cteq66"
