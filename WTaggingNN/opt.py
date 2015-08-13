def createOutYaml(config_out, job_id, filename_test):
    import os
    if not os.path.exists('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/trained/'):
        os.makedirs('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/trained/')
    outName = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/trained/config-trained-'+str(job_id)+'.yaml'
    f_out = open(outName,'w')
    f_in = open('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/template.yaml','r')
    file_holder = '  file: '+filename_test
    
    for l in f_in:
        if l.find('file:') != -1:
            f_out.write(file_holder+'\n')
        elif l.find('TAGGERYAML') != -1:
            f_out.write(l.replace('TAGGERYAML', config_out))
        else:
            f_out.write(l)
            
    f_out.close()
    f_in.close()
    return outName
    

def objective(filename_train, filename_test, tree, config, learning, momentum, regularize, uepochs, sepochs, formula, job_id):
    import subprocess
    import sys, time
    import WTaggingPerf as ww
    import yaml
    import os
    formula_out = '--formula='+formula+''
    print formula_out
    if not os.path.exists('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/output_config'):
        os.makedirs('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/output_config')
    save = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/output_config/config-output-'+str(job_id)+'.yaml'
    args = ['/Disk/ecdf-nfs-ppe/atlas/users/tibristo/AGILEPack/AGILEPackTrainer','--file='+filename_train,'--tree='+tree,'--shuffle','--save='+str(save),'--config='+config,'--learning='+str(learning),'--momentum='+str(momentum),'--regularize='+str(regularize),'--uepochs='+str(uepochs),'--sepochs='+str(sepochs),'--batch=1',formula_out]
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    lines_iter = iter(popen.stdout.readline,b"")
    f_log = open('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/db_logs/currentout.log','w')
    for line in lines_iter:
        f_log.write(line+'\n')
        f_log.flush()
        print line
    f_log.close()
    #popen.wait()

    outName = createOutYaml(save, job_id, filename_test)
    
    f = open(outName, 'r')

    schema = yaml.load(f)
    
    taggers = ww.generate_taggers(schema)
    
    auc = ww.plot_roc(taggers, schema, 'ROC_plots_' + outName.replace('.yaml','') + '.pdf', save_arr = True, logscale=False)
    # minimise this
    print auc
    return 1/auc

def main(job_id, params):
    from math import exp
    filename ='/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/'
    #filename = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/BosonTagging/'
    # these are the files used for the mva. Use one for training and one for testing.
    #persist/data_features_l_2_10__000.root
    #    persist/data_features_l_2_10__001.root
    #        persist/data_features_l_2_10__002.root

    #filename_test = filename+'persist/datatest_testroot_001.root'
    #filename_train = filename+'persist/datatrain_testroot_001.root'

    filename_test=filename+'folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_500_1000_nomw_mergedtest_cv_001.root'
    filename_train=filename+'folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedM_loose_v2_500_1000_nomw_mergedtrain_cv_001.root'
    #filename_train = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedL_ranged_v2_1000_1500_nomw_mergedtrain_cv_001.root'
    #filename_test = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedL_ranged_v2_1000_1500_nomw_mergedtest_cv_001.root'
    #config = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/tim-config.yaml'
    config = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/tim-config_full.yaml'
    formula= 'label~*-thrustmin-thrustmaj-yfilt-angularity-foxwolfram20-tau21-pt-m-eta-phi-tauwta2-tauwta1-tau2-tau1-split12 | weight'
    #formula= 'label~*-ThrustMin-ThrustMaj-YFilt-Tau21-Dip12|weight'
    auc = objective(filename_train, filename_test, 'outputTree', config, exp(params['log_learning'][0]), params['momentum'][0], exp(params['log_regularize'][0]), params['uepochs'][0], params['sepochs'][0], formula, job_id)
    print 'Job_id %d' % job_id
    print params
    print auc
    return auc


if __name__=="__main__":
    params = {'log_learning':[-2],'momentum':[0.75],'log_regularize':[-6],'uepochs':[120],'sepochs':[140]}

    main(1,params)
