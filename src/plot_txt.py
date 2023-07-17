import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np
import re

import util

def get_date_obs(obsid):
    pha2 = util.pha2_file(obsid)
    return util.read_header(pha2)['date-obs'][:10], pha2

def get_stuff(args):
    obsid, pha2, date_obs = [], [], []

    files = glob.glob(args.indir+'/[0-9]*_{}_S1.txt'.format(args.type))
    for o in [ re.search('([0-9]+)', name).group() for name in files]:
        d, p = get_date_obs(o)
        if d > '2008-12-31':
            obsid.append(o)
            date_obs.append(d)
            pha2.append(p)
    sorti = sorted(range(len(date_obs)), key=lambda i: date_obs[i])
    obsid = [obsid[i] for i in sorti]
    pha2 = [pha2[i] for i in sorti]
    date_obs = [date_obs[i] for i in sorti]
    return obsid, pha2, date_obs

def get_qes(args):
    obsid, pha2, date_obs = get_stuff(args)

    wav = { 'S1' : None, 'S2' : None, 'S3' : None }
    qe = { 'S1' : [], 'S2' : [], 'S3' : [] }

    for i in range(len(obsid)):
        for j in range(1,4):
            chip = 'S' + str(j)
            filename = '{}/{}_{}_{}.txt'.format(args.indir, obsid[i], args.type, chip)
            w, q = np.loadtxt(filename, unpack=True)

            qe[chip].append(q)
            if wav[chip] is None:
                wav[chip] = w
    return obsid, wav, qe, date_obs

def plot_qes(args, obsid, wav, qe, date_obs):

    axes = plt.figure(figsize=(8.5,11))

    ci = np.linspace(1, 0, len(obsid)-1)
    ANGSTROM, LAMBDA = "Åλ"

    for i in range(1,4):
        chip = 'S'+str(i)
        ax = plt.subplot(310+i)
        for j in range(1, len(obsid)-1):
            ax.plot(wav[chip], qe[chip][j]/qe[chip][0], color=plt.cm.RdYlBu(ci[j]), label=date_obs[j])
            ax.set_ylim(0.7, 1.1)
            ax.set_title('HRC-'+chip)
            ax.set_ylabel('ARF ratio')
            ax.set_xlabel('{} ({})'.format(LAMBDA, ANGSTROM))
            #ax.set_yscale('log')
        if i==1: ax.legend()

    plt.tight_layout()
    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser('Plot HRC-S ardlib QE for each HZ43 observation.')
    parser.add_argument('-o', '--outfile', help='Output plot filename')
    parser.add_argument('-d', '--indir', default='./misc/txt', help='Input directory filename. Default is ./misc/txt')
    parser.add_argument('type', help='Input filename substrings to match')
    args = parser.parse_args()

    obsid, wav, qe, date_obs = get_qes(args)
    plot_qes(args, obsid, wav, qe, date_obs)


if __name__ == '__main__':
    main()
