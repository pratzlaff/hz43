import numpy as np
import argparse
import math
import astropy.io.fits
import glob
import re
import util

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)
#rcParams.update({'font.size': 22})

def both_detnams(o1, o2):
    return detnams(o1), detnams(o2)

def detnams(obsids):
    detnam = None
    for obsid in obsids:
        detnam_ = util.detnam(obsid)
        #detnam_ = util.read_header(util.pha2_file(obsid))['detnam']
        if detnam is None:
            detnam = detnam_
        if detnam != detnam_:
            raise ValueError("DETNAM={} is in the obsid {} header, which does not match DETNAM={} in previous headers".format(detnam_, obsid, detnam))
    return detnam

def parse_obsids(args):

    p = re.compile('\d+')

    o1 = p.findall(args.obsids1)
    o2 = p.findall(args.obsids2)

    return o1, o2

def main():

    parser = argparse.ArgumentParser(
        description='Calculate flux ratios vs wavelength for two sets of obsids'
    )
    parser.add_argument('--combine', type=int, default=128)
    parser.add_argument('--w1', type=int, default=55)
    parser.add_argument('--w2', type=int, default=73)
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('obsids1')
    parser.add_argument('obsids2')
    args = parser.parse_args()

    o1, o2 = parse_obsids(args)
    obsids = [o1, o2]

    det1, det2 = detnams(o1), detnams(o2)

    spectra = [
        util.get_spectra(o1, 'LEG', args.combine),
        util.get_spectra(o2, 'LEG', args.combine),
    ]

    labels = { 'pos':'TG\_M = +1', 'neg':'TG\_M = -1'}

    # for order in spectra[0].keys():
    #     for i in xrange(2):
    #         print(order, i)
    #         x, y, yerr = spectra[i][order]['lambda'], spectra[i][order]['flux'], spectra[i][order]['flux_err']
    #         plt.plot(x, y, label=labels[order])
    #         plt.errorbar(x, y, yerr, ecolor='k', fmt='none')
    #         plt.show()

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize=(11, 8.5))

    w1, w2 = args.w1, args.w2

    for spectrum in spectra:
        for order in spectrum:
            mask = np.logical_and(spectrum[order]['lambda']>=w1, spectrum[order]['lambda']<=w2)
            for s in 'lambda', 'flux', 'rate', 'flux_err':
                spectrum[order][s] = spectrum[order][s][mask]
        print(spectrum['pos']['lambda'])

    for i in xrange(len(spectra)):
        for order in spectra[i]:
            for s in 'flux', 'rate', 'flux_err':
                spectra[i][order][s] /= len(obsids[i])
#            print(i, order)
#            plt.plot(spectra[i][order]['lambda'], spectra[i][order]['rate'])
#            plt.show()

    ratio = {}
    ratio_err = {}
    for order in spectra[0].keys():
        ratio[order] = spectra[0][order]['flux'] / spectra[1][order]['flux']
        ratio_err[order] = ratio[order] * np.sqrt((spectra[0][order]['flux_err'] / spectra[0][order]['flux'])**2 + (spectra[1][order]['flux_err'] / spectra[1][order]['flux'])**2)

        x, y, yerr = spectra[0][order]['lambda'], ratio[order], ratio_err[order]
        print(y)
        plt.errorbar(x, y, yerr, label=labels[order])

    plt.title(r'$\textrm{All On-Axis Observations of HZ 43 With HRC/LETG}$')
    plt.legend(loc='upper right', fontsize=12)
    plt.ylabel(r'$\textrm{Flux Ratio: ' + det1 + ' / ' + det2 + r'}$')
    plt.xlabel(r'$\textrm{Wavelength (\AA)}$')
    #plt.ylim(0.8, 1.4)
    #plt.ylim((0,2))

    plt.tight_layout()

    if args.pdf:
        pdf.savefig(fig)
        pdf.close()
    else:
        plt.show()

if __name__ == '__main__':
    main()
