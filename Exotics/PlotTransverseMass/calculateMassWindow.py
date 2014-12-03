from ROOT import *
import sys
import os
import copy
import numpy as n


def setupHistogram(fname, algorithm):
    '''
    Read in the jet mass from the input file and fill it into a histogram.  The histogram gets used to calculate the mass window.
    Keyword args:
    fname --- input file name.
    algorithm --- Algorithm being used.
    '''
    # set up a struct so that the SetBranchAddress method can be used.
    # normally one would use tree.variable to access the variable, but since
    # this variable name changes for each algorithm this is a better way to do it.
    gROOT.ProcessLine("struct jet_t { Float_t mass;} ")
    # create an object of this struct
    jet = jet_t()
    
    # open file and ttree
    f = TFile.Open(fname,"READ")
    tree = f.Get("physics")
    
    # set branch address for the groomed jet mass
    tree.SetBranchAddress("jet_"+algorithm+"_m", AddressOf(jet,'mass'))
    # histogram with 300 bins and 0 - 3 TeV range
    hist = TH1F("mass","mass",300,0,300*1000)    
    # maximum number of entries in the tree
    entries = tree.GetEntries()

    # loop through
    for e in xrange(entries):
        tree.GetEntry(e)
        # fill the hist
        hist.Fill(jet.mass)
    
    return hist

    
def Qw(histo, frac=0.68):
    '''
    Method for calculating the mass window from a histogram of masses.
    Keyword args:
    histo --- input TH1F histogram of masses
    frac --- the mass window fraction. Normally use 68%
    '''
    # set up the variables that store the best found window
    minWidth = 100000.0;
    topEdge = 0.0;
    botEdge = 0.0
    maxidx = 99
    minidx = 0

    # info on histogram - number of bins and the integral
    Nbins = histo.GetNbinsX()
    integral = histo.Integral();

    # loop through each bin of the histogram
    for i in xrange(0,Nbins):

        tempFrac = 0.0
        # want to make sure we don't change i when changing imax
        imax = copy.deepcopy(i)

        # loop through until the tempFrac is above the frac (68%) criteria,
        # but making sure not to go out of range.                                                                                                                                             
        while(tempFrac<frac and imax != Nbins):
            #fraction in bin imax=0,1,2,...                                                                                                                                                     
            tempFrac+=histo.GetBinContent(imax)/integral;
            imax+=1;

        width = histo.GetBinCenter(imax) - histo.GetBinCenter(i);
        # by applying this we say that the window we have just calculate MUST have at least 68%.
        if tempFrac >= frac and width<minWidth:
            # set up the best found mass window variables
            minWidth = width;
            topEdge = histo.GetBinCenter(imax);
            botEdge = histo.GetBinCenter(i)
            minidx = copy.deepcopy(i)
            maxidx = copy.deepcopy(imax)

    return minWidth, topEdge, botEdge, minidx, maxidx


def run(fname, algorithm):
    '''
    Method for running over a single algorithm and calculating the mass window.
    Keyword args:
    fname --- the input file name
    algorithm --- the name of the algorithm
    '''

    # setup the histogram - read in the mass entries from the input file and put them into a histogram
    hist = setupHistogram(fname, algorithm)
    # calculate the width, top and bottom edges and the indices for the 68% mass window
    wid, topedge, botedge, minidx, maxidx = Qw(hist, 0.68)

    # folder where input file is
    folder = fname[:fname.rfind('/')+1]
    # write this information out to a text file
    fout = open(folder+algorithm+"_masswindow.out",'w')
    fout.write("width: "+ str(wid)+'\n')
    fout.write("top edge: "+ str(topedge)+'\n')
    fout.write("bottom edge: "+ str(botedge)+'\n')
    fout.write("minidx: "+ str(minidx)+'\n')
    fout.write("maxidx: "+ str(maxidx)+'\n')
    fout.close()


def runAll():
    '''
    Method for running over a collection of algorithms and calculating the mass window.  
    The algorithms are stored in a text file which is read in.  Each line gives the input folder, with the folder of the form ALGORITHM_fileID/.
    '''
    # open input file
    f = open("nc29_opt_v3_files.txt","read")
    # list that will store all of the filenames and corresponding algorithms
    fnames = []
    # read in file
    for l in f:
        # remove whitespace
        l = l.strip()
        # The input folder is given on each line
        folder = l
        # now get the files in this folder and iterate through, looking for the signal file
        fileslist = os.listdir(folder)
        sigfile = ""
        for fil in fileslist:
            # found the signal file!
            if fil.endswith("sig.root"):
                sigfile = folder+fil
        # the algorithm can be obtained from the folder by removing the leading
        # directory information and removing the file identifier
        alg = l.split('/')[-2][:-len('_nc29_opt_v3')]
        # add these to the list to run over
        fnames.append([sigfile,alg])

    # loop through all of the algorithms and find the mass windows
    for a in fnames:
        print "running " + a[1]
        run(a[0],a[1])

if __name__=="__main__":
    runAll()
