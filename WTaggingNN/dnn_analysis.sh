# plotting
python create_folds.py false true --key=mc15_jz5_bkg_v2 --algorithm=AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz5_v1_800_1200_mw_merged --weight=True


# creating folds
python create_folds.py true false --key=jz5 --algorithm=AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz5_v1_800_1200_mw_merged --fulldataset=AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_jz5_notcleaned_v1_800_1200_mw_merged.csv

# run the training
python opt.py

# to test on another dataset
# can change the name of the files if needed
#for x in `ls *CHANGEME*` ; do copyfile="${x/CHANGEME/NEWNAME}" ; cp $x $copyfile ; done
cp -r trained/*mc15_jz5_bkg_v2* trained_full/
cd trained_full

# if running over another dataset and they names need to be different:
# create the training files for this -> just copy the old ones and give them a new name
for x in `ls output_config/*v2.1*` ; do copyfile="${x/v2.1/v2.4}" ; cp $x $copyfile ; done

# edit full_config to have the correct testing filename, the correct sigeff and the correct bkgeff
# the signal and bkg eff come from the dataset that is being tested.
# extra settings that have now been included in the config files:
# weight_validation = 0 or 1 -> This is if we want to weight the samples we're testing on
# transform_weight_validation 0 or 1 -> This is if we want to apply the same transformation when running the testing samples as we did when training.  If false, then just use the standard weights in the file
# need to change `weight: weight` in the config file to `weight: weight_train` for the transformed weights.
# need to change the output_config file to that above

# next edit the findbests.py file and change the file_id to match the key used above and f_out to something descriptive of the sample
# then run
python findbests.py --input-folder=trained_full/ --output-id=mc15_v2.20_bkg_v1 --file-id=mc15_v2.2_bkg_v1

# Change the following in plotEvaluation.py
# key = the key from above when creating the folds or running findbests.py
# filter_id = anything you do not want to be matched in the filename. Used for filtering
# now run plotEvaluation.py
python plotEvaluation.py --key=mc15_v2.2_bkg_v1 --filter=jz5_bkg_v1 --recreate-csv

