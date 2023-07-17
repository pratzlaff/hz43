import response
import numpy as np
import argparse
import math
import astropy.io.fits
import glob
import os
import re

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)
rcParams.update({'font.size': 14})

import util

np.seterr(divide='ignore', invalid='ignore')

def main():

    parser = argparse.ArgumentParser(
        description='Calculate HRC-I/LETG HZ 43 predicted zeroth order count rate'
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('-s', '--sobsid', default='16375', help='HRC-S obsid to use. Default is 16375')
    parser.add_argument('-i', '--iobsid', default='19847', help='HRC-I ARF to use. Default is 19847')
    parser.add_argument('-m', '--minfrac', type=float, default=1, help='HRC-S response cutoff, as a percentage of maximum. Bins below this threshold will be omitted. Default is 1.')
    parser.add_argument('-e', '--emin', type=float, default=80, help='Low energy threshold, in eV. Bins below this threshold will be omitted. Default is 80.')
    args = parser.parse_args()

    i_bin_lo, i_bin_hi, i_response = response.read_arf(response.zeroth_arf_file(args.iobsid))
    i_wav = 0.5 * (i_bin_lo + i_bin_hi)

    spectrum = util.get_spectra((args.sobsid,), 'LEG', 1)

    i_response = np.interp(spectrum['pos']['lambda'], i_wav[::-1], i_response[::-1])

    w1, w2 = util.w2e(3.0), util.w2e(args.emin/1000.)
    for order in spectrum:
        mask = (spectrum[order]['lambda']>=w1) & (spectrum[order]['lambda']<=w2)
        for s in 'lambda', 'flux', 'rate', 'flux_err', 'response':
            spectrum[order][s] = spectrum[order][s][mask]

    i_response = i_response[mask]

    resp = spectrum['pos']['response'] + spectrum['neg']['response']

    flux = (spectrum['pos']['rate'] + spectrum['neg']['rate']) / resp
    wav = spectrum['pos']['lambda']

    minfrac=args.minfrac/100.
    i = (spectrum['neg']['response']>minfrac*spectrum['neg']['response'].max()) & (spectrum['pos']['response']>minfrac*spectrum['pos']['response'].max())

    i = resp > (minfrac * resp.max())

    print('Excluding %d / %d bins below the %g%% RESPONSE threshold (%.3g cm^2)' % ((i==0).sum(), i.size, args.minfrac, minfrac * resp.max()))

    wav = wav[i] ; flux = flux[i] ; i_response = i_response[i]

    print 'Predicted HRC-I/LETG 0th order count rate = %.2f' % (flux * i_response).sum()

    if args.pdf:
        pdf = PdfPages(args.pdf)

    fig = plt.figure(figsize = (8.5, 11))

    plot_dims = (2,1)
    plt.subplot2grid(plot_dims, (0,0))

    plt.plot(wav, i_response)
    plt.xlabel(r'$\lambda$ (\AA)')
    plt.ylabel(r'Effective Area ($\textrm{cm}^2$)')
    plt.title('HRC-I/LETG 0th Order EA, Obsid {}'.format(args.iobsid))

    plt.subplot2grid(plot_dims, (1,0))

    plt.plot(wav, flux)
    plt.title('HRC-S/LETG flux, obsid {}'.format(args.sobsid))
    plt.xlabel(r'$\lambda$ (\AA)')
    plt.ylabel(r'Flux (ph/sec/$\textrm{cm}^2$/bin)')

    plt.tight_layout()

    if args.pdf:
        pdf.savefig(fig)
        pdf.close()
    else:
        plt.show()

    plt.close()

if __name__ == '__main__':
    main()
