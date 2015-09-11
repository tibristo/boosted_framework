# first change the test file to the full file
for x in `ls *paramID*.yaml` ; do sed -i '/test_cv*root/c\  file: /Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_notcleaned_v2_200_1000_mw_mergedscalefull.root' $x ; done
# change the signal eff
for x in `ls *paramID*.yaml` ; do sed -i '/signaleff/c\sigeff: 0.6839' $x ; done
# change the background eff
for x in `ls *paramID*.yaml` ; do sed -i '/bkgeff/c\bkgeff: 0.0914' $x ; done

