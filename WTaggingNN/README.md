WTaggingNN
==============
<em>Python Scripts for W Tagging Performance plotting</em>
Requirements: numpy, sklearn, pandas, root_numpy and ROOT.


Usage is simple, and can be run interactively. First, put `export PYTHONPATH+=:/path/to/WTaggingNN` in your `~/.bashrc`. 

This is based off Luke's TopTagger here found here: `git clone git@github.com:lukedeo/TopTaggingPerf.git`

There are a few different things in this package.  In order to create data splits use `create_folds.py`.  This will
take as input csv files containing the data.  This is fairly easy to create if you don't have this already.  Just
use root_numpy to write a tree to csv.  The folds get stored under the folder persist/ and have the name AlgorithmBLAH_cv_xyz.root.
The folds are created with the StratifiedKFold method from sklearn.

`create_folds.py` also has a plotting method (plotFiles) to create plots of all the variables in the cv splits and to create
stats files - the mean, std, number of signal and bkg events etc.


There is a wrapper script called `opt.py` which will run the AGILEPackTrainer binary and train a neural network.  This requires
you to have already set up the c++ version of AGILEPack. Opt.py will call this in the `objective` method with a number of arguments and wait for the training
to finish.  Once this is done, it will create a yaml file for the net just trained.  It takes a template yaml file `template.yaml` and edits it
putting in the name of the test file, the training file, the algorithm and the network training parameters.

This yaml file then loaded and then the taggers are created and plotted with WTaggingNN/ROC.py and /plot_schema.py. This will return the 1.0/(background rejection power at 50% signal efficiency).

Objective is meant to return a value which can be minimised by some calling function.  Originally it was setup so that Spearmint could continuously call it with different parameters,
using a bayesian optimisation procedure.  This is still possible - the config file for this is called `config.json` and then Spearmint will call main.

Running opt.py will call runCV() which will set up a pool of 4 worker threads.  It creates a grid of parameters to go through.  All the input folds matching the given key (set in the .py file)
are used as input.  The formula and config file for AGILEPack are set and then grid_search is called, which sets up a bunch of tasks which will be sent to the thread pool.

Some things to be aware of - there is a base config file called tim-config_full.yaml.  This has the variables in the input files (the cv splits).  Spearmint does better at
using log of the learning rate and log of the regularizer - see exp(log_learning) in main()


When processing the plots with ROC.py and plot_schema.py there are a lot of things that are done.

The yaml file containing the info from the trained neural net is read in with agilepy.client.  This is then used to create predictions.
```python
>>> import WTaggingPerf as ww
>>> import agilepy.client as agilepy
>>>
>>> net = agilepy.NeuralNet()
>>> net.load('NetworkFile.yaml')
>>> predictions = net.predict(data)['label_predicted']
>>>
>>> roc = ROC_plotter(taggers)
```
