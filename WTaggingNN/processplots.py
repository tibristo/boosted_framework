#!/usr/bin/env python
import sys
import yaml
import WTaggingPerf as tt
import datetime
import os

#if __name__ == '__main__':
folder = str(sys.argv[1])
files = [f for f in os.listdir(folder) if f.endswith('.yaml'))
filenames = str(sys.argv[1])
f = open(filename, 'r')
	
schema = yaml.load(f)

	taggers = tt.generate_taggers(schema)

	now = datetime.datetime.now()
	timestamp = now.strftime("%Y-%m-%d_%H:%M")

	tt.plot_roc(taggers, schema, 'ROC_plots_' + timestamp + '.pdf', save_arr = True, logscale=False)
	#tt.plot_roc(taggers, schema, 'ROC_plots_' + timestamp + '_zoom_.pdf', 0.4, 0.6, False)
