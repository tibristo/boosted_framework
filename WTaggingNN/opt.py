import time


results = []
def logFinished(result):
    global results
    results.append(result)


def createOutYaml(config_out, job_id, filename_test, filename_train, algorithm, params={}):
    import os
    if not os.path.exists('/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/trained/'):
        os.makedirs('/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/trained/')
    outName = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/trained/config-trained-'+str(job_id)+'.yaml'
    f_out = open(outName,'w')
    f_in = open('/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/template.yaml','r')
    file_holder = '  file: '+filename_test
    
    for l in f_in:
        if l.strip().startswith('file:'):
            f_out.write(file_holder+'\n')
        elif l.find('TAGGERYAML') != -1:
            f_out.write(l.replace('TAGGERYAML', config_out))
        elif l.find('TRAINFILE') != -1:
            f_out.write('trainfile: '+filename_train+'\n')
        elif l.find('ALGORITHM') != -1:
            f_out.write('algorithm: ' + algorithm+'\n')
            # write out the parameters here
            for k in params.keys():
                f_out.write(k+': '+str(params[k])+'\n')
        else:
            f_out.write(l)
    # store the training parameters and the train file in here as well!
    # add algorithm name
    f_out.close()
    f_in.close()
    return outName
    

def objective(filename_train, filename_test, tree, config, learning, momentum, regularize, uepochs, sepochs, formula, job_id, algorithm):
    import subprocess
    import sys, time
    import WTaggingPerf as ww
    import yaml
    import os
    formula_out = '--formula='+formula+''
    print formula_out
    if not os.path.exists('/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/output_config'):
        os.makedirs('/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/output_config')
    save = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/output_config/config-output-'+str(job_id)+'.yaml'
    args = ['/Disk/ds-sopa-group/PPE/atlas/users/tibristo/AGILEPack/AGILEPackTrainer','--file='+filename_train,'--tree='+tree,'--shuffle','--save='+str(save),'--config='+config,'--learning='+str(learning),'--momentum='+str(momentum),'--regularize='+str(regularize),'--uepochs='+str(uepochs),'--sepochs='+str(sepochs),'--batch=1',formula_out]
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    lines_iter = iter(popen.stdout.readline,b"")
    for line in lines_iter:
        pass
    #    print line
    popen.wait()
    
    # the number of epochs is not stored in the yaml weights file, so we're going to save it
    # in the output one here
    params = {'uepochs':str(uepochs), 'sepochs':str(sepochs)}
    outName = createOutYaml(save, job_id, filename_test, filename_train, algorithm, params)
    
    f = open(outName, 'r')

    schema = yaml.load(f)
    
    taggers = ww.generate_taggers(schema)

    # outName has the full path in it, so the string must be split on / first.
    outName_spl = outName.split('/')[-1]
    
    rejection_power = ww.plot_roc(taggers, schema, 'ROC_plots/ROC_plots_' + outName_spl.replace('.yaml','') + '.pdf', save_arr = True, logscale=False)
    # minimise this
    print rejection_power
    return 1.0/rejection_power


def grid_search(folder, cv_split_filenames, parameter_grid, config, algo, id_tag = 'cv', formula = ''):
    """Launch all grid search evaluation tasks."""
    from sklearn.grid_search import ParameterGrid
    from math import exp
    all_tasks = []
    all_parameters = list(ParameterGrid(parameter_grid))

    print len(all_parameters)
    for i, params in enumerate(all_parameters):
        task_for_params = []
      
        for j, cv_split_filename in enumerate(cv_split_filenames):
            filename_train = cv_split_filename[0]
            filename_test = cv_split_filename[1]
            #x = objective(filename_train, filename_test, 'outputTree', config, exp(params['log_learning']), params['momentum'], exp(params['log_regularize']), params['uepochs'], params['sepochs'], formula, job_id='paramID_'+str(i)+id_tag+'ID_'+str(j), algorithm=algo)
            # Having a strange issue when trying to run this with iPython. Going to switch to using multiprocessing.Pool
            # This is how I currently run the simple tagger. See the scanBkgRej.py file in the BosonTagging repo.
            # store the parameters that will get sent to objective
            obj_params = [filename_train, filename_test, 'outputTree', config, exp(params['log_learning']), params['momentum'], exp(params['log_regularize']), params['uepochs'], params['sepochs'], formula, 'paramID_'+str(i)+id_tag+'ID_'+str(j),algo]
            #t = lb_view.apply(
            #    objective, filename_train, filename_test, 'outputTree', config, exp(params['log_learning']), params['momentum'], exp(params['log_regularize']), params['uepochs'], params['sepochs'], formula, job_id='paramID_'+str(i)+id_tag+'ID_'+str(j), algorithm=algo)
            #raw_input()
            #task_for_params.append(t) 
            all_tasks.append(obj_params)
        #all_tasks.append(task_for_params)
        
    return all_parameters, all_tasks



def find_bests(all_parameters, all_tasks, n_top=5, save=False, bests_tag='cv'):
    """Compute the mean score of the completed tasks"""
    mean_scores = []
    param_id = 0
    for param, task_group in zip(all_parameters, all_tasks):
        scores = [t.get() for t in task_group if t.ready()]
        if len(scores) == 0:
            continue
        mean_scores.append((np.mean(scores), param, param_id))
        param_id+=1
    bests = sorted(mean_scores, reverse=True, key=lambda x: x[0])[:n_top]        
    if save:
        f = open('bests/bests'+bests_tag+'.txt','w')
        for b in bests:
            f.write('mean_score: ' + str(b[0]) + ' params: ' + str(b[1]) + ' param id: ' + str(b[2])+'\n')
        f.close()
    return bests

    
def printProgress(tasks):
    total = len(tasks)
    print total
    finished = 0.0
    for t in tasks:
        if t.ready():
            finished+=1
            t.get()
    print 'finished: ' + str(finished)
    return float(finished/total)


def main(job_id, params):
    from math import exp
    filename ='/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/'
    #filename = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/BosonTagging/'
    # these are the files used for the mva. Use one for training and one for testing.
    #persist/data_features_l_2_10__000.root
    #    persist/data_features_l_2_10__001.root
    #        persist/data_features_l_2_10__002.root

    #filename_test = filename+'persist/datatest_testroot_001.root'
    #filename_train = filename+'persist/datatrain_testroot_001.root'

    filename_test=filename+'folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_200_1000_mw_mergedtest_cv_001.root'
    filename_train=filename+'folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_200_1000_mw_mergedtrain_cv_001.root'
    #filename_train = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/nnbosontagging/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedL_ranged_v2_1000_1500_nomw_mergedtrain_cv_001.root'
    #filename_test = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/nnbosontagging/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedL_ranged_v2_1000_1500_nomw_mergedtest_cv_001.root'
    #config = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/nnbosontagging/tim-config.yaml'
    config = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/tim-config_full.yaml'
    formula= 'label~*-thrustmin-thrustmaj-yfilt-angularity-foxwolfram20-pt-m-eta-phi-mu12-zcut12-planarflow-tauwta1-tau2-tau1-sphericity-eec_d2_1-tauwta2tauwta1| weight' #-tau21-tau1-tau2 <- these are not in mc15
    #formula= 'label~*-ThrustMin-ThrustMaj-YFilt-Tau21-Dip12|weight'
    algorithm = 'AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_200_1000_mw'
    auc = objective(filename_train, filename_test, 'outputTree', config, exp(params['log_learning'][0]), params['momentum'][0], exp(params['log_regularize'][0]), params['uepochs'][0], params['sepochs'][0], formula, job_id, algorithm)
    print 'Job_id %d' % job_id
    print params
    print auc
    return auc


def runCV():

    from pprint import pprint
    from collections import OrderedDict
    import numpy as np
    import os
    from multiprocessing import Pool 
    #params = OrderedDict([
    #    ('log_learning', np.linspace(-7, -2, 5, dtype=int)),
    #    ('momentum', np.linspace(0.7, 0.85, 5)),
    #    ('log_regularize', np.linspace(-10, 0, 5, dtype=int)),
    #    ('uepochs',np.linspace(20, 60, 4, dtype=np.int)),
    #    ('sepochs',np.linspace(40, 80, 4, dtype=np.int))
    #])
    '''
    params = OrderedDict([
        ('log_learning', np.linspace(-10, -2, 5, dtype=int)),
        ('momentum', np.linspace(0.7, 0.85, 2)),
        ('log_regularize', np.linspace(-10, -7, 2)),
        ('uepochs',np.linspace(20, 60, 2, dtype=np.int)),
        ('sepochs',np.linspace(40, 80, 2, dtype=np.int))
    ])
    '''
    # params for the best 5 that we got using mc15 on jz5 background WITHOUT weighting the validation sample
    '''
    params = [{'momentum':[0.85],'log_regularize':[-10],'log_learning':[-4],'uepochs':[60],'sepochs':[80]},
              {'momentum':[0.85],'log_regularize':[-10],'log_learning':[-6],'uepochs':[20],'sepochs':[80]},
              {'momentum':[0.7],'log_regularize':[-10],'log_learning':[-4],'uepochs':[60],'sepochs':[80]},
              {'momentum':[0.7],'log_regularize':[-10],'log_learning':[-4],'uepochs':[20],'sepochs':[80]},
              {'momentum':[0.85],'log_regularize':[-10],'log_learning':[-4],'uepochs':[20],'sepochs':[80]}]
    '''
    # params for the best 5 that we got using mc15 on jz5 background WITH weighting the validation sample
    params = [{'momentum':[0.85],'log_regularize':[-7],'log_learning':[-6],'uepochs':[20],'sepochs':[80]},
              {'momentum':[0.85],'log_regularize':[-7],'log_learning':[-6],'uepochs':[60],'sepochs':[80]},
              {'momentum':[0.85],'log_regularize':[-7],'log_learning':[-6],'uepochs':[60],'sepochs':[40]},
              {'momentum':[0.85],'log_regularize':[-7],'log_learning':[-6],'uepochs':[20],'sepochs':[40]},
              {'momentum':[0.7],'log_regularize':[-7],'log_learning':[-4],'uepochs':[20],'sepochs':[40]}]
    # key to look for in the filenames
    #key = 'matchedM_loose_v2_200_1000_mw'
    key = '13tev_mc15_nTrk_v1'
    # output file id
    file_id = 'mc15_nTrk_v1_bkg_v1'
    folder ='/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/'
    # we need the files to be the output of what we would normally get back from
    # the cross_validation method from create_folds.py
    all_files = [f for f in os.listdir('folds/') if f.find(key) != -1 and f.endswith('.root')]
    
    full_dataset = ''
    for f in all_files:
        if f.find('full') != -1:
            all_files.remove(f)
            full_dataset = os.getcwd()+'/folds/' + f
            break
    # now we need to match all of the cross validation folds together.
    # get the training files
    train_files = filter(lambda x: x.find('train')!=-1, all_files)
    # cv folds
    filenames = map(lambda x: [folder+'folds/'+x, folder+'folds/'+x.replace('train','test')], train_files)
    
    print filenames
    config = '/Disk/ds-sopa-group/PPE/atlas/users/tibristo/boosted_framework/WTaggingNN/tim-config-nTrk.yaml'
    # training on mc15 with variables: aplanarity, eec_c2_1, split12, tauwta2tauwta1, eec_d2_1, zcut12, planarflow. However, zcut12 and split12 are highly
    # correlated and planarflow and sphericity are actually higher than eec_d2_1 when looking at the bdt feature importances!
    # for mc15_jz5_v1
    #formula= 'label~*-thrustmin-thrustmaj-yfilt-angularity-foxwolfram20-pt-m-eta-phi-mu12-sphericity-tauwta2-tauwta1| weight' # tau2, tau1, tau21 not in mc15
    # for mc15_jz5_v2 and v3
    #formula= 'label~*-thrustmin-thrustmaj-yfilt-angularity-foxwolfram20-pt-m-eta-phi-mu12-tauwta2-tauwta1-zcut12-weight| weight_train' # tau2, tau1, tau21 not in mc15
    # for mc15 with nTrack variable
    formula= 'label~*-thrustmin-thrustmaj-yfilt-angularity-foxwolfram20-pt-m-eta-phi-mu12-tauwta2-tauwta1-zcut12-weight| weight_train' # tau2, tau1, tau21 not in mc15
    # for mc15_blah_v4 -> v4 is NO WEIGHTS USED FOR TRAINING
    #formula= 'label~*-thrustmin-thrustmaj-yfilt-angularity-foxwolfram20-pt-m-eta-phi-mu12-tauwta2-tauwta1-zcut12-weight-weight_train' # tau2, tau1, tau21 not in mc15
    #formula= 'label~*-thrustmin-thrustmaj-yfilt-angularity-foxwolfram20-tau21-pt-m-eta-phi-tauwta2-tauwta1-tau2-tau1| weight'
    #algorithm = 'AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_v3_400_1600_mw'
    algorithm = 'AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_mc15_nTrk_v1_800_1200_mw'

    #allparms, alltasks = grid_search(
    #    lb_view, folder, filenames, params, config, algorithm, id_tag=file_id, formula=formula)
    allparms, alltasks = grid_search(folder, filenames, params, config, algorithm, id_tag=file_id, formula=formula)

    
    # now use multiprocessing.Pool to map these        
    pool = Pool(4)
    pool_results = []
    for a in alltasks:
        pool_results.append(pool.apply_async(objective, a, callback=logFinished))
    pool.close()
    
    #pool.join()
    prog = printProgress(pool_results)
    while prog < 1:
        time.sleep(10)
        prog = printProgress(pool_results)
        print "progress: " + str(prog)


if __name__=="__main__":
    #params = {'log_learning':[-2],'momentum':[0.75],'log_regularize':[-6],'uepochs':[120],'sepochs':[140]}
    runCV()
    #main(1,params)
