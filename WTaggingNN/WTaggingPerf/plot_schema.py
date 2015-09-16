import yaml
from .ROC import *
import client as apy
from rootpy.io import root_open
import root_numpy as rn
import numpy as np
import os
import hashlib
import copy

def _tree_to_array(schema, to_npy = False):
	print 'Loading File...'
	f = root_open(schema['sample']['file'])
	T = f[schema['sample']['tree']]
	if schema['sample'].has_key('selection') == False:
		this_sel = ''#None
	else:
		this_sel = schema['sample']['selection']
	if schema['sample'].has_key('step'):
		this_step = schema['sample']['step']
	else:
		this_step = ''#None
	print 'Pulling Tree...'
	arr = rn.tree2array(T, selection = this_sel, step = this_step)
	if to_npy is True:
		print 'Writing to *.npy file...'
		varlist = "".join(this_sel.split()).replace('(', '').replace(')', '').split('&&')
		varlist.sort()
		hash_name = os.path.basename(schema['sample']['file']) + schema['sample']['tree'] + ''.join(varlist) + str(this_step)
		m = hashlib.sha1()
		m.update(hash_name)
		np.save(os.path.dirname(schema['sample']['file']) + '/' + m.hexdigest() + '.npy', arr)
	print 'Done.'
	return arr

def _get_data_hash(schema):
	if schema['sample'].has_key('selection') == False:
		this_sel = ''#None
	else:
		this_sel = schema['sample']['selection']
	if schema['sample'].has_key('step'):
		this_step = schema['sample']['step']
	else:
		this_step = ''#None
	varlist = "".join(this_sel.split()).replace('(', '').replace(')', '').split('&&')
	varlist.sort()
	hash_name = os.path.basename(schema['sample']['file']) + schema['sample']['tree'] + ''.join(varlist) + str(this_step)
	m = hashlib.sha1()
	m.update(hash_name)
	return m.hexdigest()

def _get_data(schema):
	hashed_file = os.path.dirname(schema['sample']['file']) + '/' + _get_data_hash(schema) + '.npy'
	if os.path.isfile(hashed_file):
		print 'Matching Schema hash found! Loading from backup.'
		return np.load(hashed_file)
	else:
		print 'No matching Schema hash found. Loading from ROOT file.'
		return _tree_to_array(schema, True)


def generate_taggers(schema, tagger_name='tagger'):
	data = _get_data(schema)

	taggers = {}

	if schema.has_key('taggers'):
		for taggerfile, specifications in schema['taggers'].iteritems():
			print 'Working on ' + taggerfile
			net = apy.NeuralNet()
			net.load(taggerfile)
			predictions = net.predict(data)[0]
                        print net.predict(data)['label_predicted']
			opt = 'false'
                        # set the tagger_name to that given in the yaml
                        tagger_name = specifications['name']
                        tagger_variables = net.inputs
                        algorithm = 'AntiKt10LCTopoTrimmedPtFrac5SmallR20'
                        if schema.has_key('algorithm'):
                                algorithm = schema['algorithm']
                        
                        if schema.has_key('trainfile'):
                                train_file = schema['trainfile']
                        else:
                                train_file = 'folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_200_1000_mw_mergedtrain_cv_001.root'

                        # can we get the predictions on the training file as well? it's probably
                        # a good idea, especially for over/under fitting validation curves
                        # copy the schema file and just set the filename
                        schema_copy = copy.copy(schema)
                        schema_copy['sample']['file'] = train_file
                        data_train = _get_data(schema_copy)
                	predictions_train = net.predict(data_train)[0]
                        
                        params = {'learning_rate': net.learning_rate, 'momentum':net.momentum, 'regularize':net.regularize}
                        print params
                        # check if the schema has the epochs
                        if schema.has_key('uepochs'):
                                params['uepochs'] = int(schema['uepochs'])
                        if schema.has_key('sepochs'):
                                params['sepochs'] = int(schema['sepochs'])
                                
                        print tagger_variables
                        if schema.has_key('sigeff'):
                                sig_eff = float(schema['sigeff'])
                        else:
                                sig_eff=1.0
                        if schema.has_key('bkgeff'):
                                bkg_eff = float(schema['bkgeff'])
                        else:
                                bkg_eff=1.0
			if specifications.has_key('optimise'):
				opt = specifications['optimise']
			if schema.has_key('weightfiles') and schema['weightfiles'] == 'true':
				add_tagger(specifications['name'], specifications['color'], 
					   general_roc(data, predictions['label_predicted'], 100, name=tagger_name, signal_eff=sig_eff,bkg_eff=bkg_eff,variables=tagger_variables,params=params,weights=data[schema['weight']], tagger_file=taggerfile, train_file = train_file, algorithm=algorithm, data_train = data_train, discriminant_train=predictions_train['label_predicted']), taggers, opt)
			else:
				add_tagger(specifications['name'], specifications['color'], 
					   general_roc(data, predictions['label_predicted'], 100, name=tagger_name, signal_eff=sig_eff,bkg_eff=bkg_eff,variables=tagger_variables,params=params, tagger_file=taggerfile,train_file=train_file, algorithm=algorithm, data_train = data_train, discriminant_train=predictions_train['label_predicted']), taggers, opt)#10000), taggers)


	if schema.has_key('benchmarks'):
		if schema['benchmarks'].has_key('scans'):
			for var, specifications in schema['benchmarks']['scans'].iteritems():
				print 'Applying scan on ' + var
				if schema.has_key('weightfiles') and schema['weightfiles'] == 'true':
					add_tagger(specifications['name'], specifications['color'], 
						   general_roc_weighted(data, data[var], data['weight'], 10000), taggers)
				else:
					add_tagger(specifications['name'], specifications['color'], 
						   general_roc(data, data[var], 10000), taggers)
		if schema['benchmarks'].has_key('taggers'):
			for function, specifications in schema['benchmarks']['taggers'].iteritems():
				print 'Applying tagger defined by function ' + function
				call = "add_tagger(specifications['name'], specifications['color'], " + function + "(data, data['mcevt_weight'], 10000), taggers)"
				eval(call)

	return taggers

def plot_roc(dictionary, schema, name = None, min_eff = 0, max_eff = 1, logscale = True, save_arr = False):
	if name is None:
		m = hashlib.sha1()
		m.update(dictionary.__repr__())
		savename = 'ROC/ROC_' + m.hexdigest() + '.pdf'
	else:
		savename = name
	rejection_power, roc = ROC_plotter(dictionary, min_eff = min_eff, max_eff = max_eff, linewidth=2.1, signal = schema['signal'], background = schema['background'], title = schema['title'], logscale = logscale, save_arr = save_arr, inputfile=savename)
	
	roc.savefig(savename)
	return rejection_power


















	




	
