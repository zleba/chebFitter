def getDictForBDT(wantBig = True):
	BDTDict = [
                { "var" : 'Mbc',               "form" : 'Mbc',                               "ran" : (5.2, 5.3),    "log" : False, "loc_leg" : "best"}, 
		{ "var" : 'deltaE',            "form" : 'deltaE',                            "ran" : (-0.5, 0.5),   "log" : False, "loc_leg" : "best"},
		{ "var" : 'chiProb',           "form" : 'chiProb',                           "ran" : (-1, 1),       "log" : False, "loc_leg" : "best"},
		{ "var" : 'foxWolframR1',      "form" : 'foxWolframR1',                      "ran" : (0, 0.3),      "log" : False, "loc_leg" : "best"},
		{ "var" : 'foxWolframR2',      "form" : 'foxWolframR2',                      "ran" : (0, 1),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'foxWolframR3',      "form" : 'foxWolframR3',                      "ran" : (0, 1),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'foxWolframR4',      "form" : 'foxWolframR4',			     "ran" : (0, 1),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'HMT0',              "form" : 'harmonicMomentThrust0', 	     "ran" : (0.0, 2.0),    "log" : False, "loc_leg" : "best"},
		{ "var" : 'HMT1',              "form" : 'harmonicMomentThrust1', 	     "ran" : (-0.5, 0.5),   "log" : False, "loc_leg" : "best"},
		{ "var" : 'HMT2',              "form" : 'harmonicMomentThrust2',             "ran" : (0., 1.),      "log" : False, "loc_leg" : "best"},
		{ "var" : 'HMT3',              "form" : 'harmonicMomentThrust3',             "ran" : (-0.5, 0.5),   "log" : False, "loc_leg" : "best"},
		{ "var" : 'CCT1',              "form" : 'cleoConeThrust1',                   "ran" : (0, 10),       "log" : False, "loc_leg" : "best"},
		{ "var" : 'CCT2',              "form" : 'cleoConeThrust2',                   "ran" : (0, 10),       "log" : False, "loc_leg" : "best"},
		{ "var" : 'CCT3',              "form" : 'cleoConeThrust3',		     "ran" : (0, 8),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'CCT4',              "form" : 'cleoConeThrust4',                   "ran" : (0, 8),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'CCT5',              "form" : 'cleoConeThrust5',                   "ran" : (0, 6),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'CCT6',              "form" : 'cleoConeThrust6',                   "ran" : (0, 5),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'CCT7',              "form" : 'cleoConeThrust7',                   "ran" : (0, 4),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'thrustAxisCosTheta',"form" : 'thrustAxisCosTheta',                "ran" : (-1, 1),       "log" : False, "loc_leg" : "best"},
		{ "var" : 'cosTBTO',           "form" : 'cosTBTO',                           "ran" : (0., 1.),      "log" : False, "loc_leg" : "best"},
		{ "var" : 'cosTBz',            "form" : 'cosTBz',                            "ran" : (0., 1.1),     "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW02',            "form" : "KSFWVariables__bohso02__bc",        "ran" : (-1., 1.),     "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW04',            "form" : "KSFWVariables__bohso04__bc",        "ran" : (-1., 1.),     "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW10',            "form" : "KSFWVariables__bohso10__bc",        "ran" : (0., 1.5),     "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW12',            "form" : "KSFWVariables__bohso12__bc",        "ran" : (-0.5, 1.),    "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW14',            "form" : "KSFWVariables__bohso12__bc",        "ran" : (-0.5, 1.),    "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW20',            "form" : "KSFWVariables__bohso20__bc",        "ran" : (0., 1.),      "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW22',            "form" : "KSFWVariables__bohso22__bc",        "ran" : (-1, 1.),      "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW24',            "form" : "KSFWVariables__bohso24__bc",        "ran" : (-1, 1.),      "log" : False, "loc_leg" : "best"}]

	extraDict = [
                { "var" : 'HMT4',              "form" : 'harmonicMomentThrust4',             "ran" : (-1, 1),       "log" : False, "loc_leg" : "best"},
		{ "var" : 'CCT8',              "form" : 'cleoConeThrust8',                   "ran" : (0, 2),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW00',            "form" : "KSFWVariables__bohso00__bc",        "ran" : (0, 3.),       "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW0',             "form" : "KSFWVariables__bohoo0__bc",         "ran" : (0., 1.),      "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW1',             "form" : "KSFWVariables__bohoo1__bc",         "ran" : (-0.1, 0.1),   "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW2',             "form" : "KSFWVariables__bohoo2__bc",         "ran" : (-0.05, 0.2),  "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW3',             "form" : "KSFWVariables__bohoo3__bc",         "ran" : (-0.05, 0.05), "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW4',             "form" : "KSFWVariables__bohoo4__bc",         "ran" : (-0.05, 0.1),  "log" : False, "loc_leg" : "best"},
		{ "var" : 'sphericity',        "form" : 'sphericity',                        "ran" : (0, 1),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'aplanarity',        "form" : 'aplanarity',                        "ran" : (0, 0.3),      "log" : False, "loc_leg" : "best"},
		{ "var" : 'thrust',            "form" : 'thrust',                            "ran" : (0.5, 1.0),    "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW_bomm2',        "form" : "KSFWVariables__boet__bc",           "ran" : (-20, 20.),    "log" : False, "loc_leg" : "best"},
		{ "var" : 'KSFW_et',           "form" : "KSFWVariables__boet__bc",           "ran" : (0., 15.),     "log" : False, "loc_leg" : "best"}]

	candidates = [
		{ "var" : 'NCandidates',       "form" : '__ncandidates__',                   "ran" : (0, 5),        "log" : False, "loc_leg" : "best"},
		{ "var" : 'candidate',         "form" : '__candidate__',                     "ran" : (0, 5),        "log" : True,  "loc_leg" : "best"},
                { "var" : 'experiment',        "form" : '__experiment__',                    "ran" : (0, 10000),    "log" : False, "loc_leg" : "best"},
                { "var" : 'run',               "form" : '__run__',                           "ran" : (0, 10000),    "log" : False, "loc_leg" : "best"},
                { "var" : 'event',             "form" : '__event__',                         "ran" : (0, 10000000), "log" : False, "loc_leg" : "best"},
                { "var" : 'chiProbRank',       "form" : 'chiProb_rank',                      "ran" : (0, 10),       "log" : False, "loc_leg" : "best"},
                { "var" : 'randomRank',        "form" : 'random_rank',                       "ran" : (0, 10),       "log" : True,  "loc_leg" : "best"},
                { "var" : 'absDeltaERank',     "form" : 'absDeltaE_rank',                    "ran" : (0, 10),       "log" : True,  "loc_leg" : "best"},
                { "var" : 'cosMdstIndRank',    "form" : 'cosMdstInd_rank',                   "ran" : (0, 10),       "log" : True,  "loc_leg" : "best"},
                { "var" : 'dilutionRank',      "form" : 'dilution_rank',                     "ran" : (0, 10),       "log" : True,  "loc_leg" : "best"}]

	if wantBig: return BDTDict + extraDict + candidates

	return BDTDict


def daugterNames(channel):
	if channel == "B0toPiDPtoK2Pi":   return ['Kp1FD', 'pim1FD', 'pim2FD', 'prongFB']
	if channel == "B0toPiDStoK2Pi":   return ['Kp1FDFDst', 'pim1FDFDst', 'pimFDst', 'prongFB']
	if channel == "B0toPiDStoK4Pi":   return ['Kp1FDFDst', 'pim1FDFDst', 'pim2FDFDst', 'pip1FDFDst', 'pimFDst', 'prongFB']
	if channel == "BPtoPiD0toKPi":    return ['Kp1FD', 'pim1FD', 'prongFB']
	if channel == "BPtoPiD0toK3Pi":   return ['Kp1FD', 'pim1FD', 'pim2FD', 'pip1FD', 'prongFB']
	if channel == "BPtoJPsiKtoEE":    return ['lepp', 'lepm', 'Kp']
	if channel == "BPtoJPsiKtoMuMu":  return ['lepp', 'lepm', 'Kp']
	if channel == "B0toJPsiKStoMuMu": return ['lepp', 'lepm', 'pip1FK', 'pim1FK']
	if channel == "B0toJPsiKStoEE":   return ['lepp', 'lepm', 'pip1FK', 'pim1FK']

def compositeNames(channel):
        if channel == "B0toPiDPtoK2Pi":   return ['D']
        if channel == "B0toPiDStoK2Pi":   return ['Dst', 'DFDst']
        if channel == "B0toPiDStoK4Pi":   return ['Dst', 'DFDst']
        if channel == "BPtoPiD0toKPi":    return ['D']
        if channel == "BPtoPiD0toK3Pi":   return ['D']
        if channel == "BPtoJPsiKtoEE":    return ['JPsi']
        if channel == "BPtoJPsiKtoMuMu":  return ['JPsi']
        if channel == "B0toJPsiKStoMuMu": return ['JPsi', 'K_S0']
        if channel == "B0toJPsiKStoEE":   return ['JPsi', 'K_S0']

def getDictForBToDPi(channel):

	#prepare the variables for the daughter tracks
	daughterDicts = [ 
		{"var" : "_kaonID",              "form" : "", "ran" : (0, 1.),           "log" : False, "loc_leg" : "best"},
	        {"var" : "_pionID",              "form" : "", "ran" : (0., 1.),          "log" : False, "loc_leg" : "best"}, 
	        {"var" : "_electronID",          "form" : "", "ran" : (0., 1.),          "log" : False, "loc_leg" : "best"},
                {"var" : "_muonID",              "form" : "", "ran" : (0., 1.),          "log" : False, "loc_leg" : "best"},
                {"var" : "_nMatchedKLMClusters", "form" : "", "ran" : (0., 5.),          "log" : False, "loc_leg" : "best"},
                {"var" : "_clusterE",            "form" : "", "ran" : (0., 3.),          "log" : False, "loc_leg" : "best"},
                {"var" : "_clusterEoP",          "form" : "", "ran" : (0., 2.),          "log" : False, "loc_leg" : "best"},
                {"var" : "_nPXDHits",            "form" : "", "ran" : (0., 4.),          "log" : False, "loc_leg" : "best"},
                {"var" : "_nSVDHits",            "form" : "", "ran" : (0., 20.),         "log" : False, "loc_leg" : "best"},
                {"var" : "_nCDCHits",            "form" : "", "ran" : (0.,100.),         "log" : False, "loc_leg" : "best"},
                {"var" : "_cms_p",               "form" : "", "ran" : (0., 4.),          "log" : False, "loc_leg" : "best"},
		{"var" : "_cms_pt",              "form" : "", "ran" : (0., 4.),          "log" : False, "loc_leg" : "best"},
                {"var" : "_p",                   "form" : "", "ran" : (0., 5.),          "log" : False, "loc_leg" : "best"},
	        {"var" : "_pt",                  "form" : "", "ran" : (0., 4.),          "log" : False, "loc_leg" : "best"},
	        {"var" : "_pValue",              "form" : "", "ran" : (0, 1),            "log" : False, "loc_leg" : "best"},
	        {"var" : "_theta",               "form" : "", "ran" : (0, 3.0),          "log" : False, "loc_leg" : "best"},
		{"var" : "_phi",                 "form" : "", "ran" : (-3.1416, 3.1416), "log" : False, "loc_leg" : "best"},
		{"var" : "_charge",              "form" : "", "ran" : (-1.5, 1.5),       "log" : False, "loc_leg" : "best"}]

	totDicts = []

	if channel == "B0toJPsiKStoMuMu" or channel == "B0toJPsiKStoEE":
		for daughterName in daugterNames(channel):
			for dictio in daughterDicts:
				totDicts.append({"var" : daughterName+dictio["var"], "form": daughterName+dictio["var"], "ran" : dictio["ran"], "log" : dictio["log"], "loc_leg" : dictio["loc_leg"]})
		totDicts.append({"var" : "pip1FK_d0", "form" : "pip1FK_d0", "ran" : (-10., 10.), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "pip1FK_z0", "form" : "pip1FK_z0", "ran" : (-10., 10.), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "pip1FK_dr", "form" : "pip1FK_dr", "ran" : (0.00, 10.), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "pip1FK_dz", "form" : "pip1FK_dz", "ran" : (-10., 10.), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "pim1FK_d0", "form" : "pim1FK_d0", "ran" : (-10., 10.), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "pim1FK_z0", "form" : "pim1FK_z0", "ran" : (-10., 10.), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "pim1FK_dr", "form" : "pim1FK_dr", "ran" : (0.00, 10.), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "pim1FK_dz", "form" : "pim1FK_dz", "ran" : (-10., 10.), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "lepp_d0",   "form" : "lepp_d0", "ran" : (-0.07, 0.07),     "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "lepp_z0",   "form" : "lepp_z0", "ran" : (-0.2, 0.2),       "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "lepp_dr",   "form" : "lepp_dr", "ran" : (0, 0.02),         "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "lepp_dz",   "form" : "lepp_dz", "ran" : (-0.01, 0.01),     "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "lepm_d0",   "form" : "lepm_d0", "ran" : (-0.07, 0.07),     "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "lepm_z0",   "form" : "lepm_z0", "ran" : (-0.2, 0.2),       "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "lepm_dr",   "form" : "lepm_dr", "ran" : (0, 0.02),         "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "lepm_dz",   "form" : "lepm_dz", "ran" : (-0.01, 0.01),     "log" : False, "loc_leg" : "best"})

	else: 
		daughterDicts.append({"var" : "_d0",                  "form" : "", "ran" : (-0.07, 0.07),     "log" : False, "loc_leg" : "best"})
		daughterDicts.append({"var" : "_z0",                  "form" : "", "ran" : (-0.2, 0.2),       "log" : False, "loc_leg" : "best"})
		daughterDicts.append({"var" : "_dr",                  "form" : "", "ran" : (0, 0.02),         "log" : False, "loc_leg" : "best"})
		daughterDicts.append({"var" : "_dz",                  "form" : "", "ran" : (-0.01, 0.01),     "log" : False, "loc_leg" : "best"})
		for daughterName in daugterNames(channel):
			for dictio in daughterDicts:
				totDicts.append({"var" : daughterName+dictio["var"], "form": daughterName+dictio["var"], "ran" : dictio["ran"], "log" : dictio["log"], "loc_leg" : dictio["loc_leg"]})

	#prepare the variables for composite particles
	compDicts = [
		{"var" : "_cms_p",   "form" : "", "ran" : (0., 4.),    "log" : False, "loc_leg" : "best"},
		{"var" : "_cms_pt",  "form" : "", "ran" : (0., 4.),    "log" : False, "loc_leg" : "best"},
		{"var" : "_p",       "form" : "", "ran" : (0., 5.),    "log" : False, "loc_leg" : "best"},
		{"var" : "_pt",      "form" : "", "ran" : (0., 4.),    "log" : False, "loc_leg" : "best"},
		{"var" : "_chiProb", "form" : "", "ran" : (0., 1.),    "log" : False, "loc_leg" : "best"},
		{"var" : "_charge",  "form" : "", "ran" : (-1.5, 1.5), "log" : False, "loc_leg" : "best"}]

	if channel == "B0toJPsiKStoMuMu" or channel == "B0toJPsiKStoEE":
		for compositeName in compositeNames(channel):
                	for dictio in compDicts:
                        	totDicts.append({"var" : compositeName+dictio["var"], "form": compositeName+dictio["var"], "ran" : dictio["ran"], "log" : dictio["log"], "loc_leg" : dictio["loc_leg"]})

		totDicts.append({"var" : "K_S0_x",                      "form" : "K_S0_x",                      "ran" : (-20.0, 20.0),  "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_y",                      "form" : "K_S0_y",                      "ran" : (-20.0, 20.0),  "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_z",                      "form" : "K_S0_z",                      "ran" : (-20.0, 20.0),  "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_x_uncertainty",          "form" : "K_S0_x_uncertainty",          "ran" : (0.0, 0.2),     "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_y_uncertainty",          "form" : "K_S0_y_uncertainty",          "ran" : (0.0, 0.2),     "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_z_uncertainty",          "form" : "K_S0_z_uncertainty",          "ran" : (0.0, 0.4),     "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_flightTime",             "form" : "K_S0_flightTime",             "ran" : (0., .5),       "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_flightDistance",         "form" : "K_S0_flightDistance",         "ran" : (0., 20.),      "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_distance",               "form" : "K_S0_distance",               "ran" : (0., 20.),      "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_significanceOfDistance", "form" : "K_S0_significanceOfDistance", "ran" : (0., 20.),      "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_dx",                     "form" : "K_S0_dx",                     "ran" : (-20., 20.),    "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_dy",                     "form" : "K_S0_dy",                     "ran" : (-20., 20.),    "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_dz",                     "form" : "K_S0_dz",                     "ran" : (-20., 20.),    "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_dr",                     "form" : "K_S0_dr",                     "ran" : (0., 10.),      "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_dphi",                   "form" : "K_S0_dphi",                   "ran" : (-3.141,3.141), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_dcosTheta",              "form" : "K_S0_dcosTheta",              "ran" : (-1.,1.),       "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_difX",                   "form" : "K_S0_x - IPX",                "ran" : (-20., 20.),     "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_difY",                   "form" : "K_S0_y - IPY",                "ran" : (-20., 20.),     "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "K_S0_difZ",                   "form" : "K_S0_z - IPZ",                "ran" : (-20., 20.),     "log" : False, "loc_leg" : "best"})

		totDicts.append({"var" : "JPsi_x",                      "form" : "JPsi_x",                      "ran" : (-.1, .1),      "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_y",                      "form" : "JPsi_y",                      "ran" : (-.1, .1),      "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_z",                      "form" : "JPsi_z",                      "ran" : (-.2, .2),      "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_x_uncertainty",          "form" : "JPsi_x_uncertainty",          "ran" : (0.0, 0.005),   "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_y_uncertainty",          "form" : "JPsi_y_uncertainty",          "ran" : (0.0, 0.005),   "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_flightTime",             "form" : "JPsi_flightTime",             "ran" : (-0.005,0.005), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_flightDistance",         "form" : "JPsi_flightDistance",         "ran" : (-0.1,0.1),     "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_distance",               "form" : "JPsi_distance",               "ran" : (0.0,0.2),      "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_significanceOfDistance", "form" : "JPsi_significanceOfDistance", "ran" : (0.0,5.0),      "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_dx",                     "form" : "JPsi_dx",                     "ran" : (-0.03,0.03),   "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_dy",                     "form" : "JPsi_dy",                     "ran" : (-0.03,0.03),   "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_dz",                     "form" : "JPsi_dz",                     "ran" : (-0.1,0.1),     "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_dr",                     "form" : "JPsi_dr",                     "ran" : (0.,0.03),      "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_dphi",                   "form" : "JPsi_dphi",                   "ran" : (-3.141,3.141), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "JPsi_dcosTheta",              "form" : "JPsi_dcosTheta",              "ran" : (-1.,1.),       "log" : False, "loc_leg" : "best"})
	
	else:
		compDicts.append({"var" : "_x",                      "form" : "", "ran" : (-.1, .1),      "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_y",                      "form" : "", "ran" : (-.1, .1),      "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_z",                      "form" : "", "ran" : (-.2, .2),      "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_x_uncertainty",          "form" : "", "ran" : (0., 0.005),    "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_y_uncertainty",          "form" : "", "ran" : (0., 0.005),    "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_z_uncertainty",          "form" : "", "ran" : (0., 0.005),    "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_flightTime",             "form" : "", "ran" : (-0.005,0.005), "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_flightDistance",         "form" : "", "ran" : (-0.1, 0.1),    "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_distance",               "form" : "", "ran" : (0.0,0.2),      "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_significanceOfDistance", "form" : "", "ran" : (0.0,5.0),      "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_dx",                     "form" : "", "ran" : (-0.03,0.03),   "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_dy",                     "form" : "", "ran" : (-0.03,0.03),   "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_dz",                     "form" : "", "ran" : (-0.1,0.1),     "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_dr",                     "form" : "", "ran" : (0.,0.03),      "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_dphi",                   "form" : "", "ran" : (-3.141,3.141), "log" : False, "loc_leg" : "best"})
		compDicts.append({"var" : "_dcosTheta",              "form" : "", "ran" : (-1.,1.),       "log" : False, "loc_leg" : "best"})

		for compositeName in compositeNames(channel):
			for dictio in compDicts:
				totDicts.append({"var" : compositeName+dictio["var"], "form": compositeName+dictio["var"], "ran" : dictio["ran"], "log" : dictio["log"], "loc_leg" : dictio["loc_leg"]})

	#{"var" : "_Q",              "form" : "",                                     "ran" : (0., 0.02),   "log" : False, "loc_leg" : "best"},
	if 'D'     in compositeNames(channel): totDicts.append({"var" : "D_M",     "form" : "D_M",     "ran" : (1.80, 1.90), "log" : False, "loc_leg" : "best"})
	if 'DFDst' in compositeNames(channel): 
		totDicts.append({"var" : "DFDst_M",      "form" : "DFDst_M",         "ran" : (1.80, 1.90), "log" : False, "loc_leg" : "best"})
		totDicts.append({"var" : "dif_DstDFDst", "form" : "Dst_M - DFDst_M", "ran" : (0.14, 0.16), "log" : False, "loc_leg" : "best"})
	if 'Dst'   in compositeNames(channel): totDicts.append({"var" : "Dst_M",   "form" : "Dst_M",   "ran" : (1.95, 2.15), "log" : False, "loc_leg" : "best"})
	if 'JPsi'  in compositeNames(channel): totDicts.append({"var" : "JPsi_M",  "form" : "JPsi_M",  "ran" : (2.80, 3.20), "log" : False, "loc_leg" : "upper left"})
	if 'K_S0'  in compositeNames(channel): totDicts.append({"var" : "K_S0_M",  "form" : "K_S0_M",  "ran" : (0.45, 0.55), "log" : False, "loc_leg" : "best"})

	tagDicts = [
		{"var" : "x",                           "form" : "x",                                     "ran" : (-.1, .1),     "log" : False, "loc_leg" : "best"},
                {"var" : "y",                           "form" : "y",                                     "ran" : (-.1, .1),     "log" : False, "loc_leg" : "best"},
                {"var" : "z",                           "form" : "z",                                     "ran" : (-.2, .2),     "log" : False, "loc_leg" : "best"},
                {"var" : "x_uncertainty",               "form" : "x_uncertainty",                         "ran" : (0., 0.008),   "log" : False, "loc_leg" : "best"},
                {"var" : "y_uncertainty",               "form" : "y_uncertainty",                         "ran" : (0., 0.008),   "log" : False, "loc_leg" : "best"},
                {"var" : "z_uncertainty",               "form" : "z_uncertainty",                         "ran" : (0., 0.008),   "log" : False, "loc_leg" : "best"},
		{"var" : "charge",                      "form" : "charge",                                "ran" : (-1.5, 1.5),   "log" : False, "loc_leg" : "best"},
                {"var" : "log_dist",                    "form" : "TMath::Log10(distance)",                "ran" : (-3., 1.),     "log" : False, "loc_leg" : "best"},
                {"var" : "log_signDist",                "form" : "TMath::Log10(significanceOfDistance)",  "ran" : (-1., 1.),     "log" : False, "loc_leg" : "best"},
		{"var" : "FANN_qrCombined",             "form" : "FANN_qrCombined",                       "ran" : (-1., 1.),     "log" : False, "loc_leg" : "best"},
	        {"var" : "FBDT_qrCombined",             "form" : "FBDT_qrCombined",                       "ran" : (-1., 1.),     "log" : False, "loc_leg" : "best"},
             	{"var" : "DeltaT",                      "form" : "DeltaT",                                "ran" : (-10., 10.),   "log" : False, "loc_leg" : "best"},
                {"var" : "DeltaTErr",                   "form" : "DeltaTErr",                             "ran" : (0., 6.),      "log" : False, "loc_leg" : "best"},
             	{"var" : "DeltaZ",                      "form" : "DeltaZ",                                "ran" : (-0.1, 0.1),   "log" : False, "loc_leg" : "best"},
		{"var" : "DeltaX",                      "form" : "x - TagVx",                             "ran" : (-0.03, 0.03), "log" : False, "loc_leg" : "best"},
		{"var" : "DeltaY",                      "form" : "y - TagVy",                             "ran" : (-0.03, 0.03), "log" : False, "loc_leg" : "best"},
		{"var" : "TagVx",                       "form" : "TagVx",                                 "ran" : (-0.10, 0.10), "log" : False, "loc_leg" : "best"},
		{"var" : "TagVy",                       "form" : "TagVy",                                 "ran" : (-0.10, 0.10), "log" : False, "loc_leg" : "best"},
             	{"var" : "TagVz",                       "form" : "TagVz",                                 "ran" : (-0.15, 0.15), "log" : False, "loc_leg" : "best"},
             	{"var" : "TagVxErr",                    "form" : "TagVxErr",                              "ran" : (0, 0.01),     "log" : False, "loc_leg" : "best"},
             	{"var" : "TagVyErr",                    "form" : "TagVyErr",                              "ran" : (0, 0.01),     "log" : False, "loc_leg" : "best"},
             	{"var" : "TagVzErr",                    "form" : "TagVzErr",                              "ran" : (0, 0.01),     "log" : False, "loc_leg" : "best"},
                {"var" : "TagVNTracks",                 "form" : "TagVNTracks",                           "ran" : (0, 12),       "log" : False, "loc_leg" : "best"},
                {"var" : "TagVNFitTracks",              "form" : "TagVNFitTracks",                        "ran" : (0, 12),       "log" : False, "loc_leg" : "best"},
                {"var" : "TagVType",                    "form" : "TagVType",                              "ran" : (0, 5),        "log" : False, "loc_leg" : "best"},
		{"var" : "TagVpVal",                    "form" : "TagVpVal",                              "ran" : (0, 1),        "log" : False, "loc_leg" : "best"},
                {"var" : "TagVDistanceToConstraint",    "form" : "TagVDistanceToConstraint",              "ran" : (0, 0.1),      "log" : False, "loc_leg" : "best"},
                {"var" : "TagVDistanceToConstraintErr", "form" : "TagVDistanceToConstraintErr",           "ran" : (0, .025),     "log" : False, "loc_leg" : "best"},
		{"var" : "TagVChi2",                    "form" : "TagVChi2",                              "ran" : (.0, 30.0),    "log" : False, "loc_leg" : "best"},
		{"var" : 'nTracks',                     "form" : "nTracks",                               "ran" : (0., 30.),     "log" : False, "loc_leg" : "best"},
		# {"var" : 'nECLClusters',                "form" : "nECLClusters",                          "ran" : (0., 100.),   "log" : False, "loc_leg" : "best"},
		{"var" : 'totalPhotonECL',              "form" : "totalPhotonEECL",                       "ran" : (0., 5.),      "log" : False, "loc_leg" : "best"}]

	vexDicts = [
		{"var" : "distance",                    "form" : "distance",                              "ran" : (0.0,0.2),      "log" : False, "loc_leg" : "best"},
		{"var" : "significanceOfDistance",      "form" : "significanceOfDistance",                "ran" : (0.0,5.0),      "log" : False, "loc_leg" : "best"},
                {"var" : "dx",                          "form" : "dx",                                    "ran" : (-0.03,0.03),   "log" : False, "loc_leg" : "best"},
                {"var" : "dy",                          "form" : "dy",                                    "ran" : (-0.03,0.03),   "log" : False, "loc_leg" : "best"},
                {"var" : "dz",                          "form" : "dz",                                    "ran" : (-0.1,0.1),     "log" : False, "loc_leg" : "best"},
                {"var" : "dr",                          "form" : "dr",                                    "ran" : (0.,0.04),      "log" : False, "loc_leg" : "best"},
                {"var" : "dphi",                        "form" : "dphi",                                  "ran" : (-3.141,3.141), "log" : False, "loc_leg" : "best"},
                {"var" : "dcosTheta",                   "form" : "dcosTheta",                             "ran" : (-1.,1.),       "log" : False, "loc_leg" : "best"},
                {"var" : "difX",                        "form" : "x - IPX",                               "ran" : (-0.03,0.03),   "log" : False, "loc_leg" : "best"},
                {"var" : "difY",                        "form" : "y - IPY",                               "ran" : (-0.03,0.03),   "log" : False, "loc_leg" : "best"},
                {"var" : "difZ",                        "form" : "z - IPZ",                               "ran" : (-0.1,0.1),     "log" : False, "loc_leg" : "best"},
                {"var" : "IPX",                         "form" : "IPX",                                   "ran" : (-.1, .1),      "log" : False, "loc_leg" : "best"},
                {"var" : "IPY",                         "form" : "IPY",                                   "ran" : (-.1, .1),      "log" : False, "loc_leg" : "best"},
                {"var" : "IPZ",                         "form" : "IPZ",                                   "ran" : (-.2, .2),      "log" : False, "loc_leg" : "best"},
                {"var" : "IPXErr",                      "form" : "IPCov__bo0__cm0__bc",                   "ran" : (0., 0.000012), "log" : False, "loc_leg" : "best"},
                {"var" : "IPYErr",                      "form" : "IPCov__bo1__cm1__bc",                   "ran" : (0., 0.00002),  "log" : False, "loc_leg" : "best"},
                {"var" : "IPZErr",                      "form" : "IPCov__bo2__cm2__bc",                   "ran" : (0., 0.002),    "log" : False, "loc_leg" : "best"},
                {"var" : "IPCov01",                     "form" : "IPCov__bo0__cm1__bc",                   "ran" : (0., 0.0001),   "log" : False, "loc_leg" : "best"},
                {"var" : "IPCov02",                     "form" : "IPCov__bo0__cm2__bc",                   "ran" : (0., 0.0001),   "log" : False, "loc_leg" : "best"},
                {"var" : "IPCov12",                     "form" : "IPCov__bo1__cm2__bc",                   "ran" : (0., 0.0001),   "log" : False, "loc_leg" : "best"}]

	totDicts = totDicts + tagDicts + vexDicts

	return totDicts

def makeBDPiPlot(dictio, pdData, pdBSig, pdBBkg, pdContinuum, plt, dirName, NBins = 40, wantSignal = True, channel = "B0ToDmPip"):
	"""
	Three possible channels:
	B0ToDmPip
	BpToD0Pip2Prong
	BpToD0Pip4Prong
	"""
	import numpy as np
	import histoPandaTools as hPd
	#first, scale to correct area 
	wSigToCont = 0.02036
	if channel == "BpToD0Pip2Prong": wSigToCont = 80/9407.34
	if channel == "BpToD0Pip4Prong": wSigToCont = (1./4.) * (3.87/4.196) 
	nData = len(pdData)
	nMC = len(pdContinuum) + len(pdBBkg) + wSigToCont * len(pdBSig)
	wMC = nData/nMC
	wSig = wMC*wSigToCont
	dataHist = hPd.getSymmetricErrorHist(pdData[dictio['var']], bins=NBins, range=dictio['ran'])#, weights = [1.]*len(pdData))
	MCHist = hPd.getSymmetricErrorHist(pdBSig[dictio['var']], bins=NBins, range=dictio['ran'], 
	                                   weights = np.array([wSig]*len(pdBSig)))
	maxSig =  MCHist[1].max()
	fig, axes = plt.subplots(figsize=(10,6))
	hBkg = plt.hist([pdContinuum[dictio['var']], pdBBkg[dictio['var']]], bins=NBins, range=dictio['ran'], 
	                weights=[[wMC]*len(pdContinuum), [wMC]*len(pdBBkg)], stacked = True, 
	                label = ["$q\overline{q}$ MC", "$b\overline{b}$ MC"])
	maxBkg = hBkg[0][1].max()
	if wantSignal == True:
		hSig = plt.hist( pdBSig[dictio['var']], bins=NBins, range=dictio['ran'], 
	        	         weights = (0.75*maxBkg/maxSig)*np.array([wSig]*len(pdBSig)),
	                	 histtype = "step", lw=3, color= "red", label = "signal MC (scaled)")

	hTot = plt.hist([pdContinuum[dictio['var']], pdBBkg[dictio['var']], pdBSig[dictio['var']]], bins=NBins, 
	                range=dictio['ran'], weights=[[wMC]*len(pdContinuum), [wMC]*len(pdBBkg), [wSig]*len(pdBSig)], stacked = True, 
	               histtype = "step", color=["none", "none", "black"], label=["", "", "sig MC"])
	plt.errorbar(x=dataHist[0], y= dataHist[1], xerr=dataHist[2], yerr=dataHist[3], color="black", fmt="o", 
	             label = "data Exp7+8")
	axes.set_xlabel(dictio['var'], fontsize = 20)
	axes.set_ylabel("#events", fontsize = 20)
	axes.set_xlim(dictio['ran'])
	axes.tick_params(axis='both', which='major', labelsize=16)
	axes.legend(loc=dictio['loc_leg'], prop={'size':20})
        #axes.legend(loc="upper left", prop={'size':20})
	if dictio['log']:
	    axes.semilogy()
	plt.tick_params(labelsize = 20)
	fig.tight_layout()
	if dirName != "":
	    fig.savefig(dirName+"/"+dictio['var']+".png", bbox_inches='tight')
