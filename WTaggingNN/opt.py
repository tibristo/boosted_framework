def createOutYaml(config_out, job_id):
    outName = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/config-trained-'+str(job_id)+'.yaml'
    f_out = open(outName,'w')
    f_in = open('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/template.yaml','r')
    for l in f_in:
        if l.find('TAGGERYAML') == -1:
            f_out.write(l)
        else:
            f_out.write(l.replace('TAGGERYAML', config_out))
    f_out.close()
    f_in.close()
    return outName
    

def objective(filename, tree, config, learning, momentum, regularize, uepochs, sepochs, formula, job_id):
    import subprocess
    import sys, time
    import WTaggingPerf as ww
    import yaml
    formula_out = '--formula='+formula+''
    print formula_out
    save = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/boosted_framework/WTaggingNN/config-output-'+str(job_id)+'.yaml'
    args = ['/Disk/ecdf-nfs-ppe/atlas/users/tibristo/AGILEPack/AGILEPackTrainer','--file='+filename,'--tree='+tree,'--shuffle','--save='+str(save),'--config='+config,'--learning='+str(learning),'--momentum='+str(momentum),'--regularize='+str(regularize),'--uepochs='+str(uepochs),'--sepochs='+str(sepochs),'--batch=1',formula_out]
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    lines_iter = iter(popen.stdout.readline,b"")
    f_log = open('/Disk/ecdf-nfs-ppe/atlas/users/tibristo/db_logs/currentout.log','w')
    for line in lines_iter:
        f_log.write(line+'\n')
        f_log.flush()
        print line
    f_log.close()
    #popen.wait()

    outName = createOutYaml(save, job_id)
    
    f = open(outName, 'r')

    schema = yaml.load(f)
    
    taggers = ww.generate_taggers(schema)
    
    auc = ww.plot_roc(taggers, schema, 'ROC_plots_' + outName.replace('.yaml','') + '.pdf', save_arr = True, logscale=False)
    # minimise this
    print auc
    return 1/auc

def main(job_id, params):
    from math import exp
    filename = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/folds/AntiKt10LCTopoTrimmedPtFrac5SmallR20_13tev_matchedL_ranged_v2_1000_1500_nomw_mergedtrain_cv_001.root'
    config = '/Disk/ecdf-nfs-ppe/atlas/users/tibristo/nnbosontagging/tim-config.yaml'
    formula= 'label~*-thrustmin-thrustmaj-yfilt-mu12-angularity-aplanarity-planarflow-split12-foxwolfram20-sphericity-zcut12-eta-phi-tau21-tauwta2tauwta1-pt | weight'
    auc = objective(filename, 'outputTree', config, exp(params['log_learning'][0]), params['momentum'][0], exp(params['log_regularize'][0]), params['uepochs'][0], params['sepochs'][0], formula, job_id)
    print 'Job_id %d' % job_id
    print params
    print auc
    return auc


if __name__=="__main__":
    params = {'log_learning':[-10],'momentum':[0.75],'log_regularize':[-10],'uepochs':[20],'sepochs':[40]}

    main(1,params)
