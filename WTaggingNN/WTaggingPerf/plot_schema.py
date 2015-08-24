import yaml
from .ROC import *
import client as apy
from rootpy.io import root_open
import root_numpy as rn
import numpy as np
import os
import hashlib

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
                        if specifications.has_key('algorithm'):
                                algorithm = specifications['algorithm']
                        
                        if specifications.has_key('trainfile'):
                                train_file = specifications['trainfile']
                        else:
                                train_file = 'folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_200_1000_mw_mergedtrain_cv_001.root'
                        params = {'learning_rate': net.learning_rate, 'momentum':net.momentum, 'regularize':net.regularize}
                        print tagger_variables
                        if specifications.has_key('signaleff'):
                                sig_eff = float(specifications['signaleff'].strip())
                        else:
                                sig_eff=1.0
                        if specifications.has_key('bkgeff'):
                                bkg_eff = float(specifications['bkgeff'].strip())
                        else:
                                bkg_eff=1.0
			if specifications.has_key('optimise'):
				opt = specifications['optimise']
			if schema.has_key('weightfiles') and schema['weightfiles'] == 'true':
				add_tagger(specifications['name'], specifications['color'], 
					   general_roc(data, predictions['label_predicted'], 100, name=tagger_name, signal_eff=sig_eff,bkg_eff=bkg_eff,variables=tagger_variables,params=params,weights=data[schema['weight']], tagger_file=taggerfile, train_file = train_file, algorithm=algorithm), taggers, opt)
			else:
				add_tagger(specifications['name'], specifications['color'], 
					   general_roc(data, predictions['label_predicted'], 100, name=tagger_name, signal_eff=sig_eff,bkg_eff=bkg_eff,variables=tagger_variables,params=params, tagger_file=taggerfile,train_file=train_file, algorithm=algorithm), taggers, opt)#10000), taggers)


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


















	




	
