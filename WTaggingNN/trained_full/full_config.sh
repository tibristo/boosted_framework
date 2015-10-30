# first change the test file to the full file
#for x in `ls *bdt*.yaml` ; do sed -i '/test_cv.*root/c\  file: /Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_notcleaned_v2_200_1000_mw_mergedscalefull.root' $x ; done

#for x in `ls *jz5_v2_bkg_v1*.yaml` ; do l=`grep cv $x | head -1`; cv=${l:(-6):1}; sed -i '/test_cv.*root/c\  file: /Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz5_notcleaned_v2_800_1200_mw_merged_full_scaled_cv_00'${cv}'full.root' $x ; done

#for x in `ls *bdt*.yaml` ; do sed -i 's/notcleaned/loose/g' $x ; done
# change the signal eff
for x in `ls *jz5_v2_bkg_v1*.yaml` ; do sed -i '/sigeff/c\sigeff: 0.704' $x ; done
#for x in `ls *bdt*.yaml` ; do sed -i '/signaleff/c\sigeff: 0.6839' $x ; done
# change the background eff
#for x in `ls *bdt*.yaml` ; do sed -i '/bkgeff/c\bkgeff: 0.0914' $x ; done
for x in `ls *jz5_v2_bkg_v1*.yaml` ; do sed -i '/bkgeff/c\bkgeff: 0.1195' $x ; done

