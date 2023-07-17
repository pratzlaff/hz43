import os, errno, sys
import cPickle
import carmcmc as cm
import numpy as np
import argparse
import glob
import re

npaths=1000
nsamples=20000

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def basename(file):
    return re.compile('([^/]+)\.\w+$').findall(file)[0]

def get_samples(ratio_file):

    global nsamples

    pickle_dir = './pickle'
    pickle_file = '%s/%s_samples.bin' % (pickle_dir, basename(ratio_file))

    try:
        f = open(pickle_file, 'r')
        samples = cPickle.load(f)
        f.close()

    except:

        mkdir_p(pickle_dir)
        wav, ratio, err = np.loadtxt(ratio_file, usecols=(0,2,3), unpack=True)
        models = cm.CarmaModel(wav, ratio, err, p=2, q=0)

#        models.choose_order(10, njobs=-1) # takes 30 minutes to run on legs
        samples = models.run_mcmc(nsamples)
        f = open(pickle_file, 'w')
        cPickle.dump(samples, f, 2)
        f.close()

    return samples

def get_sims(ratio_file):

    global npaths

    samples = get_samples(ratio_file)

    year, order = re.compile('(\d{4})_([a-z]{3})').findall(ratio_file)[0]
    
    chipy_file = { 'neg' : '../chipy_3.txt', 'pos' : '../chipy_1.txt' }[order]
    chipy, wav = np.loadtxt(chipy_file, unpack=True)
    wav = np.sort(wav)

    sims_shape = ( npaths, wav.size )
    pickle_dir = './pickle'
    pickle_file = '%s/%s_sims.bin' % (pickle_dir, basename(ratio_file))

    try:
        f = open(pickle_file, 'r')
        sims = cPickle.load(f)
        f.close()

        if cmp( sims.shape, sims_shape ):
            raise "The simulations must be recalculated"

    except:
        mkdir_p(pickle_dir)
        sims = np.zeros(sims_shape)
        for i in range(npaths):
            sims[i] = samples.simulate(wav, bestfit='random')

        f = open(pickle_file, 'w')
        cPickle.dump(sims, f, 2)
        f.close()

    return wav, sims

def main():

    parser = argparse.ArgumentParser(
        description='Perform CARMA simulation models HRC-S HZ 43 model ratios vs wavelength.'
    )

    args = parser.parse_args()

    ratio_files = glob.glob('ratios/[0-9]*_*.ratio')

    for file in ratio_files:
        wav, sims = get_sims(file)
        sims_mean = sims.mean(axis=0)
        ratio_file = "%s.model" % basename(file)
        np.savetxt(ratio_file, np.transpose([wav, sims.mean(axis=0)]), fmt=['%.2f', '%f'], delimiter="\t", header="lambda\tflux_ratio")

if __name__ == '__main__':
    main()
