import numpy as np
import agilepy as apy
net = apy.NeuralNetwork()
net.load('../nnbosontagging/tim-output-noweight.yaml')
T = apy.root.tree_reader('../nnbosontagging/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedL_ranged_v2_1000_1500_nomw_mergedtest_cv_000.root')
T.get_tree('outputTree')
T.get_array(branches = net.branches)
X = T.to_array()
predictions = net.predict(X)
