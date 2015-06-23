./AGILEPackTrainer --file=/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedL_ranged_v2_1000_1500_nomw_mergedtrain_cv_001.root --tree=outputTree --shuffle --save=/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/tim-output-weight_v5.yaml --config=/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/tim-config.yaml --learning=0.0009 --momentum=0.85 --regularize=0.000001 --uepochs=10 --sepochs=40 --batch=1 --formula="label~*-thrustmin-thrustmaj-yfilt-mu12-angularity-aplanarity-planarflow-split12-foxwolfram20-sphericity-zcut12-eta-phi-tau2-tau1-pt | weight"
#./AGILEPackTrainer --file=/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedL_ranged_v2_1000_1500_nomw_mergedtrain_cv_001.root --tree=outputTree --shuffle --save=/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/tim-output-weight_v5.yaml --config=/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/tim-config.yaml --learning=0.001 --momentum=0.75 --regularize=0.001 --uepochs=40 --sepochs=80 --batch=1 --formula="label~*-thrustmin-thrustmaj-yfilt-mu12-angularity-aplanarity-planarflow-split12-foxwolfram20-sphericity-zcut12-eta-phi-tau2-tau1-pt | weight"
#--formula="label~*-ThrustMin-ThrustMax-YFilt-Mu12-Angularity-Aplanarity-PlanarFlow-SPLIT12-FoxWolfram20 | weight" 
#--start: 0 --end: 11000000 
#--struct=27 30 25 18 7 3
#  thrustmin: double
#  tau1: double
#  sphericity: double
#  foxwolfram20: double
#  tau21: double
#  thrustmaj: double
#  eec_c2_1: double
#  pt: double
#  eec_c2_2: double
#  phi: double
#  split12: double
#  tauwta2tauwta1: double
#  eec_d2_1: double
#  yfilt: double
#  mu12: double
#  tauwta2: double
#  angularity: double
#  zcut12: double
#  tau2: double
#  eec_d2_2: double
#  eta: double
#  tauwta1: double
#  planarflow: double
#  label: int
#  weight: double