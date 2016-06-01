# first change the test file to the full file



for x in `ls *mc15_nTrk_v10_bkg_v1*.yaml` ; do l=`grep cv $x | head -1`; cv=${l:(-6):1}; sed -i '/test_cv.*root/c\  file: /Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_nTrk_notcleaned_v10_400_1200_mw_merged_full_scaled_cv_00'${cv}'full.root' $x ; done
for x in `ls *mc15_nTrk_v3_bkg_v1*.yaml` ; do l=`grep cv $x | head -1`; cv=${l:(-6):1}; sed -i '/test_cv.*root/c\  file: /Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz5_nTrk_notcleaned_v3_800_1200_mw_merged_full_scaled_cv_00'${cv}'full.root' $x ; done
#for x in `ls *mc15_jz5_nTrk_v1_bkg_v3*.yaml` ; do l=`grep cv $x | head -1`; cv=${l:(-6):1}; sed -i '/test_cv.*root/c\  file: /Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz5_nTrk_notcleaned_v1_800_1200_mw_merged_full_scaled_cv_00'${cv}'full.root' $x ; done

#for x in `ls *bdt*.yaml` ; do sed -i 's/notcleaned/loose/g' $x ; done

# change the signal eff
for x in `ls *mc15_nTrk_v10_bkg_v1*.yaml` ; do sed -i '/sigeff/c\sigeff: 0.6907' $x ; done
for x in `ls *mc15_jz5_nTrk_v3_bkg_v1*.yaml` ; do sed -i '/sigeff/c\sigeff: 0.704' $x ; done
#for x in `ls *mc15_jz5_nTrk_v1_bkg_v4*.yaml` ; do sed -i '/sigeff/c\sigeff: 0.704' $x ; done
#for x in `ls *mc15_jz5_nTrk_v1_bkg_v3*.yaml` ; do sed -i '/sigeff/c\sigeff: 0.704' $x ; done


# change the background eff
for x in `ls *mc15_nTrk_v10_bkg_v1*.yaml` ; do sed -i '/bkgeff/c\bkgeff: 0.0936' $x ; done
for x in `ls *mc15_jz5_nTrk_v3_bkg_v1*.yaml` ; do sed -i '/bkgeff/c\bkgeff: 0.1195' $x ; done
#for x in `ls *mc15_jz5_nTrk_v1_bkg_v4*.yaml` ; do sed -i '/bkgeff/c\bkgeff: 0.1195' $x ; done
#for x in `ls *mc15_jz5_nTrk_v1_bkg_v3*.yaml` ; do sed -i '/bkgeff/c\bkgeff: 0.1195' $x ; done

#for x in `ls *mc15_jz5_nTrk_v1_bkg_v4*.yaml` ; do sed -i '/weightfiles:/c\weightfiles: False' $x ; done

# add in the weight-validation  and transform-weight-validation options.  Could always just stick it in on the 10th line as a quick fix, or, what's probably safer, is to just put it in after weight
#for x in `ls *mc15_v3.3_bkg_v1*.yaml` ; do sed -i '/weight-validation: True/c\weight-validation: False' $x ; done
#for x in `ls *mc15_v3.3_bkg_v1*.yaml` ; do sed -i '/transform-weight-validation: True/c\transform-weight-validation: False' $x ; done
#for x in `ls *mc15_v2.5_bkg_v1*.yaml` ; do sed -i '/weight: weight/c\weight: weight_train\nweight-validation: True' $x ; done #\ntransform-weight-validation: True' $x ; done

#for x in `ls *mc15_v2.4_bkg_v1*.yaml` ; do sed -i '/weight: weight/c\weight: weight_train' $x ; done

#for x in `ls *mc15_v3.3_bkg_v1*.yaml` ; do sed -i 's/v3.1/v3.3/g' $x ; done
#for x in `ls *mc15_jz5_nTrk_v1*.yaml` ; do sed -i 's/__cv/_cv/g' $x ; done
