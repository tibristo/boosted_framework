import numpy as np
import math
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.mlab import griddata
from matplotlib import rc
import numpy.ma as ma
from sklearn.metrics import roc_curve, auc, f1_score, precision_score, recall_score
from ROOT import TH1F, TH2D, TCanvas, TFile, TNamed, gROOT
from root_numpy import fill_hist
import functions as fn
import modelEvaluation as me
gROOT.SetBatch(True)
import os
import sys

def data_load(filename):
    data = np.genfromtxt(filename, delimiter = ',', names = True)
    return data

def pt_plot(data, pt_label = 'pt', weighted = None, log = False, bins = 100, logcount = False):
	fig = plt.figure(figsize=(11.69, 8.27), dpi=100) 
	ax = plt.subplot(1,1,1)
	if log:
		ax.set_xscale('log')
	pt = (data[pt_label]) / 1000
	bins = np.linspace(np.min(pt), np.max(pt), bins)
	if log:
		plt.xlabel(r"$\log(p_T)$ in GeV", fontsize=20)
	else:
		plt.xlabel(r"$p_T$ in GeV", fontsize=20)
	plt.ylabel(r"Count", fontsize=20)
	plt.title(r"Distribution of Jet $p_T$")
	plt.xlim(np.min(pt), np.max(pt))
	if weighted is None:
		plt.hist(pt[data['label'] == 1], histtype='step', log = logcount, bins = bins, color = 'red', label = r"$W'$ jets")
		plt.hist(pt[data['label'] == 0], histtype='step', log = logcount, bins = bins, color = 'blue', label = "JZxW Jets")
	else:
		plt.hist(pt[data['label'] == 1], histtype='step', log = logcount, bins = bins, color = 'red', label = r"$W'$ jets", weights=weighted[data['label'] == 1])
		plt.hist(pt[data['label'] == 0], histtype='step', log = logcount, bins = bins, color = 'blue', label = "JZxW Jets", weights=weighted[data['label'] == 0])
	
	plt.legend(loc = 1)
	plt.show()
	return fig


def sb_plot(data, signal = 'label', label = 'pt', disc = None, weighted = None, log = False, bins = 100, logcount = False):
	fig = plt.figure(figsize=(11.69, 8.27), dpi=100) 
	ax = plt.subplot(1,1,1)
	if disc is None:
		pt = data[label]
	else:
		pt = disc
	if log:
		ax.set_xscale('log')
	bins = np.linspace(np.min(pt), np.max(pt), bins)
	if log:
		plt.xlabel(r"$\log$"+ str(label), fontsize=20)
	else:
		plt.xlabel(r"" + str(label), fontsize=20)

	plt.ylabel(r"Count", fontsize=20)
	plt.title(r"Distribution of Jet " + label)
	
	plt.xlim(np.min(pt), np.max(pt))
	if weighted is None:
		plt.hist(pt[data[signal] == 1], histtype='step', log = logcount, bins = bins, color = 'red', label = r"Signal = " + str(signal))
		plt.hist(pt[data[signal] == 0], histtype='step', log = logcount, bins = bins, color = 'blue', label = "Background")
	else:
		plt.hist(pt[data[signal] == 1], histtype='step', log = logcount, bins = bins, color = 'red', label = r"Signal = " + str(signal), weights=weighted[data[signal] == 1])
		plt.hist(pt[data[signal] == 0], histtype='step', log = logcount, bins = bins, color = 'blue', label = "Background", weights=weighted[data[signal] == 0])
	
	plt.legend(loc = 1)
	plt.show()
	return fig


def pre_process(data):
	#data = data[data['fjet_pt']/1000 < 1100]
	#data = data[data['fjet_pt']/1000 > 550]
	#data = data[np.abs(data['fjet_eta']) < 1.2 ]
	#data = data[data['fjet_Tau2'] > -10]
	#data = data[data['fjet_Tau3'] > -10]
	return data

def score(top_ind, qcd_ind, discriminant, sample_weight = None):
    # the score of the tagger we will define as the % correct classification when taking
    # a cut at 50% signal efficiency
    # Accuracy is how many correctly classified signal AND background
    # what is the middle element's value in the signal array? that is the cut

    sig_med = np.median(discriminant[top_ind])
    bkg_med = np.median(discriminant[qcd_ind])
    top_disc = discriminant[top_ind]
    qcd_disc = discriminant[qcd_ind]
    # at 50% signal eff
    cut = top_disc[len(top_disc)/2]

    sig_correct = top_disc >= cut
    sig_incorrect = top_disc < cut
    bkg_correct = qcd_disc < cut
    bkg_incorrect = qcd_disc >= cut
    # if signal median is less than bkg, swap signal and background
    if sig_med < bkg_med:
        sig_correct, sig_incorrect = sig_incorrect, sig_correct
        bkg_correct, bkg_incorrect = bkg_incorrect, bkg_correct

    if sample_weight is None:
        sig_score = sig_correct.sum()
        sig_failed = sig_incorrect.sum()
        bkg_score = bkg_correct.sum()
        bkg_failed = bkg_incorrect.sum()
        score = float(sig_score+bkg_score)/float(discriminant.shape[0])
    else:
        sig_weights = sample_weight[top_ind]
        sig_score = np.average(sig_correct, weights=sig_weights)
        sig_failed = np.average(sig_incorrect, weights=sig_weights)
        bkg_weights = sample_weight[qcd_ind]
        bkg_score = np.average(bkg_correct, weights=bkg_weights)
        bkg_failed = np.average(bkg_incorrect, weights=bkg_weights)
        all_correct = np.concatenate((sig_correct,bkg_correct),axis=0)
        all_correct_weights = np.concatenate((sig_weights, bkg_weights),axis=0)
        score = np.average(all_correct, weights=all_correct_weights)
    # calculate the different scoring metrics
    # accuracy = correct/total -> doesn't discrim between signal and bkg
    # precision = correct signal / (TP+FP)
    # recall = TP/(TP+FN)
    # f1 = 2*precion.recall/ (precision+recall)
    # TP = 1 when 1, FP = 1 when 0, TN = 0 when 0, FN = 0 when 1
    # we need to decide what our decision value is! what probability do we cut on?
    # choose 50% since that is what we normally go for in the cut-based tagger
    # for the bkg rejection power (at 50% signal)
    accuracy = score
    # Usage for precision, recall and f1 asks for blah_score(y_true, y_pred, sample_weight=weight)
    # top_ind is effectively y_true, the ground truth target labels -> True for label==1, and False for label==0. This is set in
    # general_roc() below.
    # y_pred is an array of the estimated targets
    if sig_med < bkg_med:
        y_pred = discriminant < cut
    else:
        y_pred = discriminant >= cut

    precision = precision_score(top_ind, y_pred, sample_weight=sample_weight)

    recall = recall_score(top_ind, y_pred, sample_weight=sample_weight)
    
    f1 = f1_score(top_ind, y_pred, sample_weight=sample_weight)
    
    return cut, accuracy, precision, recall, f1
    
def general_roc(data, discriminant, bins = 2000, inverse=False, name="", signal_eff=1.0, bkg_eff=1.0, variables = [], params=[], weights_list=[], tagger_file='', train_file = '', algorithm = '', data_train = [], discriminant_train = [], weight_validation = False, tx_weight_validation = False):
	top = data[:]['label']

        train_avail = len(data_train) > 0
        
	# qcd_total = np.sum(top == 0)
	# top_total = np.sum(top == 1)
	bincount = bins

	top_ind = data[:]['label'] == 1
	qcd_ind = data[:]['label'] == 0

        if tx_weight_validation and 'weight_train' in data[:].dtype.names:
            weights = data[:]['weight_train'] if weight_validation else None # if we want to test with the transformed weights, then this needs to be "weight_train"
        elif tx_weight_validation and weight_validation:
            weights = data[:]['weight']# if weight_validation else None 
            # now transform the weights
            # all signal events have weight 1, all bkg events have weight 1/arctan(weight)
            print 'transforming weights'
            for idx in xrange(0, weights.shape[0]):
                if data[idx]['label'] == 1:
                    weights[idx] = 1.0
                else:
                    weights[idx] = 1.0/np.arctan(weights[idx])
            print 'done transforming weights'
        else:
            weights = data[:]['weight'] if weight_validation else None 
        print 'weight_validation ' + str(weight_validation)
        print 'tx_weight_validation ' + str(tx_weight_validation)
        print weights
        discriminant_bins = np.linspace(np.min(discriminant), np.max(discriminant), bins)
        if train_avail:
            top_ind_tr = data_train[:]['label'] == 1
	    qcd_ind_tr = data_train[:]['label'] == 0
            weights_train = data_train[:]['weight_train'] if weight_validation else None

        fpr,tpr,thresholds = roc_curve(top,discriminant,sample_weight=weights)

        # get the scores for the validation
        cut, accuracy, precision, recall, f1 = score(top_ind, qcd_ind, discriminant, sample_weight = weights)
        # get the scores for the training
        if train_avail:
            cut_tr, accuracy_tr, precision_tr, recall_tr, f1_tr = score(top_ind_tr, qcd_ind_tr, discriminant_train, sample_weight = weights_train)
        taggers = variables
        # model could store the net, but it's easier to just say agile and then store the
        # yaml filename
        # job_id comes from the tagger_file name
        if tagger_file != "":
            job_id = tagger_file[tagger_file.rfind('/')+1:].replace('.yaml','')
        else:
            job_id = ""
        
        model = me.modelEvaluation(fpr, tpr, thresholds, tagger_file, params, job_id, taggers, algorithm, accuracy, train_file)

        # set the scores
        model.setScores('test',accuracy=accuracy, precision=precision, recall=recall, f1=f1, cut=cut)
        if train_avail:
            model.setScores('train',accuracy=accuracy_tr, precision=precision_tr, recall=recall_tr, f1=f1_tr, cut=cut_tr)
            model.setNumberTrainEvents(signal=len(top_ind_tr),background=len(qcd_ind_tr))
            
        # get the efficiencies
        model.setSigEff(signal_eff)
        model.setBkgEff(bkg_eff)
        model.setOutputPath('ROC_root')
        model.setOutputPrefix('AGILE_')
        model.setProbas(discriminant, top_ind, qcd_ind, weights)
        model.toROOT()
        print job_id
        # now pickle the file
        pickle_filename = 'evaluationObjects/'+job_id+'.pickle'

        import pickle
        try:
            with open(pickle_filename,'w') as d:
                pickle.dump(model,d)
            d.close()
        except:
            msg = 'unable to dump ' + job_id+ ' object\n'
            msg += str(sys.exc_info()[0])
            with open(pickle_filename,'w') as d:
                pickle.dump(msg, d)
            d.close()
            print msg

        return model.ROC_sig_efficiency, model.ROC_bkg_rejection, model.hist_sig, model.hist_bkg, model.roc_graph




def ROC_plotter(taggerdict, min_eff = 0, max_eff = 1, linewidth = 1.4, pp = False, signal = "$W'", background = "JZ3-7", title = "W Tagging Efficiency vs. Rejection", logscale = True, save_arr = False, inputfile = ''):
	fig = plt.figure(figsize=(11.69, 8.27), dpi=100)
	ax = fig.add_subplot(111)
	plt.xlim(min_eff,max_eff)
	plt.grid(b = True, which = 'minor')
	plt.grid(b = True, which = 'major')
        plt.clf()
	max_ = 0
        roc_auc = 0
        rejpow = -1
        tagger_roc_auc = -1
        tagger_rejpow = -1

	for tagger, data in taggerdict.iteritems():
            
	    #sel = (data['efficiency'] >= min_eff) & (data['efficiency'] <= max_eff)
            # if using sel, then add [sel] to the end of all data[whatever] to get data[whatever][sel]
	    if np.max(data['rejection']) > max_:
		max_ = np.max(data['rejection'])
	    if save_arr:
		ar = np.zeros((data['rejection'].shape[0], 2))
		ar[:, 0] = data['efficiency']
		ar[:, 1] = data['rejection']
		np.savetxt(tagger+'_save.csv', ar, delimiter=',')

	    plt.plot(data['efficiency'], data['rejection'], '-', label = r''+tagger, color = data['color'], linewidth=linewidth)
            
            # in order to have some consistency between TaggerTim.py, mva_tools.py and this
            # we are going to use the fn.GetBGRej50() method from functions.py
            rej = fn.GetBGRej50(data['roc_curve'])
            if rej != 1:
                rejpow = 1/(1-rej)
            else:
                rejpow = -1

            print 'rejection power: ' + str(rejpow)
            tagger_rejpow = rejpow


	ax = plt.subplot(1,1,1)
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
		item.set_fontsize(20)
	
	if logscale == True:	
		plt.ylim(1,10 ** 3)
		ax.set_yscale('log')
	ax.set_xlabel(r'W efficiency')
	ax.set_ylabel(r"QCD rejection (1-eff)")

	plt.legend()
	plt.title(r''+title + ': ' + str(rejpow))
	if pp:
		pp.savefig(fig)
                return tagger_rejpow
	else:
		plt.show()
		return tagger_rejpow, fig


def add_tagger(name, color, tagger_pair, dictref, opt=False):
    if len(tagger_pair) == 2:
	dictref.update({name : {'efficiency' : tagger_pair[0], 'rejection' : tagger_pair[1], 'color' : color, 'optimise': opt, 'hist_sig':None, 'hist_bkg': None}})
    else:
        dictref.update({name : {'efficiency' : tagger_pair[0], 'rejection' : tagger_pair[1], 'color' : color, 'optimise': opt, 'hist_sig':tagger_pair[2], 'hist_bkg': tagger_pair[3], 'roc_curve':tagger_pair[4]}})








