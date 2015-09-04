import WTaggingPerf as ww
import yaml
f = open('trained_full/config-trained-10.yaml')
schema = yaml.load(f)
taggers = ww.generate_taggers(schema)
auc = ww.plot_roc(taggers, schema, 'ROC_plots_testing.pdf', save_arr = True, logscale=False)
print auc
