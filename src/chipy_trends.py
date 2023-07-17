import hz43
import response
import numpy as np
import astropy.io.fits
from astropy.io.fits.column import _parse_tdim
import argparse

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)
#rcParams.update({'font.size': 14})

import flux as flux_module
import util

maxorder=None

def get_evt2_cols(obsid, tg_reprocess='tg_reprocess'):
     hdulist = astropy.io.fits.open(util.evt2_file(obsid, tg_reprocess=tg_reprocess))
     d = hdulist['events'].data
#     return d['chip_id'], d['chipy'], d['tg_lam']
     i = np.where(np.abs(d['tg_d']) < 0.1)
     return d['chip_id'][i], d['chipy'][i], d['tg_lam'][i]

def calc_flux(src, bg, resp, bin_lo, bin_hi, hdr, factor=1):

    bg_mult = util.bg_mult(hdr)

    net = src - bg * bg_mult
    net_var = src + bg * bg_mult**2

    rate = net / hdr['exposure']
    rate_err = np.sqrt(net_var) / hdr['exposure']

    # correct for higher-order contributions, to get just 1st order
    rate *= factor
    rate_err *= factor

    flux, flux_err = flux_module.rate2flux(rate, rate_err, resp, bin_lo, bin_hi)

    return rate, rate_err, flux, flux_err

def get_spectra(obsid, resp, factor, tg_reprocess='tg_reprocess'):
    d, h = util.read_pha2(util.pha2_file(obsid, tg_reprocess=tg_reprocess))
    orders = { 'neg' : -1, 'pos' : +1 }
    spectra = { }
    
    for o in orders:
        ind = np.where(d['tg_m']==orders[o])[0][0]
        src = d['counts'][ind]
        bg = d['background_up'][ind] + d['background_down'][ind]
        rate, rate_err, flux, flux_err = calc_flux(src, bg, resp[o], d['bin_lo'][ind], d['bin_hi'][ind], h, factor[o])

        spectra[o] = {
             'bin_lo' : d['bin_lo'][ind],
             'bin_hi' : d['bin_hi'][ind],
            'rate' : rate,
            'rate_err' : rate_err,
            'flux' : flux,
            'flux_err' : flux_err
        }

    return spectra

def get_wav_ranges(obsid, args, tg_reprocess='tg_reprocess'):

    chip_id, chipy, tg_lam = get_evt2_cols(obsid, tg_reprocess=tg_reprocess)

    chipy_ranges = read_qeu_chipy(args.qeu)

    chip_ids = { 'neg' : 3, 'pos' : 1 }

    wav_ranges = {}

    for o in chipy_ranges:

        wav_ranges[o] = np.zeros(chipy_ranges[o].shape)

        ii = np.where(chip_id == chip_ids[o])[0]
        chipy_ = chipy[ii]
        tg_lam_ = tg_lam[ii]

        for i in range(chipy_ranges[o].shape[1]):
            chipy_min, chipy_max = chipy_ranges[o][0][i], chipy_ranges[o][1][i]
            jj = np.where(np.logical_and(chipy_ >= chipy_min, chipy_ < chipy_max))[0]
            if jj is not None and jj.size>90:
                wav_ranges[o][0][i] = tg_lam_[jj].min()
                wav_ranges[o][1][i] = tg_lam_[jj].max()
            else:
                wav_ranges[o][0][i] = np.nan
                wav_ranges[o][1][i] = np.nan

    return wav_ranges

def get_fracs(obsid, args, tg_reprocess='tg_reprocess'):
    global maxorder

    bin_lo, bin_hi, resp = response.get_response(obsid, 'LEG', maxorder=maxorder)
    model_wav, model_flux = hz43.model()
    wav = 0.5*(bin_lo+bin_hi)
    energy = util.w2e(wav)
    model_flux = np.interp(wav, model_wav, model_flux) * (bin_hi - bin_lo)

    # ratio of 1st-order counts / all orders
    predicted_rate = { order : model_flux / 1.602e-9 / energy * resp[order] for order in resp }
    factor = { order : predicted_rate[order][0] / predicted_rate[order].sum(axis=0) for order in resp }

    spectra = get_spectra(obsid, resp, factor, tg_reprocess='tg_reprocess')

    wav_ranges = get_wav_ranges(obsid, args)

    fracs = { o : np.zeros(wav_ranges[o].shape[1]) for o in wav_ranges }
    frac_errs = { o : np.zeros(wav_ranges[o].shape[1]) for o in wav_ranges }

    for o in fracs:

         if False and np.where(spectra[o]['bin_lo'] != bin_lo[0])[0] is not None:
              print(spectra[o]['bin_lo'].shape, bin_lo[0].shape)
              print(spectra[o]['bin_lo'], bin_lo[0])
              print(spectra[o]['bin_hi'], bin_hi[0])
              raise RuntimeError("whoops: {}".format(obsid))

         for j in range(wav_ranges[o].shape[1]):
              wav_min = wav_ranges[o][0][j]
              wav_max = wav_ranges[o][1][j]
              ind = np.where(np.logical_and(bin_lo[0]>=wav_min, bin_hi[0]<=wav_max))
              pred = predicted_rate[o][0][ind].sum()
              obs = spectra[o]['rate'][ind].sum()
              obs_err = np.sqrt((spectra[o]['rate_err'][ind]**2).sum())
              fracs[o][j] = obs/pred
              frac_errs[o][j] = obs_err/pred

    return wav_ranges, fracs, frac_errs

def calc_chipy_ranges(tdim3, _3crpx3, _3crvl3, _3cdlt3):
    nchipy = _parse_tdim(tdim3)[0]
    chipy = (np.arange(nchipy, dtype=np.float) + (_3crpx3-1)) * _3cdlt3 + _3crvl3

    ranges = np.zeros((2, nchipy), dtype=np.float)
    ranges[0] = chipy - _3cdlt3/2
    ranges[1] = chipy + _3cdlt3/2

    return ranges

def read_qeu_chipy(qeu):
    hdulist = astropy.io.fits.open(qeu)

    chipy = []

    for i in range(1,4):
        d = hdulist[i].data

        # FIXME: this picks up 101, 201, 301, which is fine for S1 and
        # S3, but not S2. Luckily, we don't care
        ii = np.where(d['regionid'] == i*100+1)
        j = ii[0][0]
        chipy.append(
            calc_chipy_ranges(d['tdim3'][j],
                         d['3crpx3'][j],
                         d['3crvl3'][j],
                         d['3cdlt3'][j]
            )
        )
    return { 'pos' : chipy[0], 'neg' : chipy[2] }
                     
def main():
     global maxorder

     parser = argparse.ArgumentParser(
          description='HRC-S/LETG HZ 43 observations deviations from model, for the RAWY regions in the QEU',
     )
     parser.add_argument('-q', '--qeu', default='/data/legs/rpete/flight/qeu/ARD/v11/hrcsD1999-07-22qeuN0011pre.fits', help='QEU file')
     parser.add_argument('-p', '--pdf', help='Save plot to named file.')
     parser.add_argument('-m', '--maxorder', help='Maximum ARF/RMF order to read.', default=3, type=int)
    
     args = parser.parse_args()

     maxorder = args.maxorder

     obsids, dates = hz43.obsids_years('HRC-S')

     chipy_ranges = read_qeu_chipy(args.qeu)

     fracs = { o : np.zeros((obsids.size, chipy_ranges[o].shape[1])) for o in chipy_ranges }
     frac_errs = { o : np.zeros((obsids.size, chipy_ranges[o].shape[1])) for o in chipy_ranges }
     wav_ranges = { o : np.zeros((obsids.size, 2, chipy_ranges[o].shape[1])) for o in chipy_ranges }

     for i in range(obsids.size):
          wav_range, frac, frac_err = get_fracs(obsids[i], args)
          for o in frac:
               fracs[o][i] = frac[o]
               frac_errs[o][i] = frac_err[o]
               wav_ranges[o][i] = wav_range[o]

     if args.pdf:
          pdf = PdfPages(args.pdf)
          fig = plt.figure(figsize = (8.5, 11))

     if False:
          plot_dims = (3, 2)

          for i in range(obsids.size):

               row = int(i/plot_dims[1]) % plot_dims[0]
               col = i % plot_dims[1]
               plt.subplot2grid(plot_dims, (row,col))

               for o in fracs:
                    x = 0.5*(chipy_ranges[o][0]+chipy_ranges[o][1])
                    y = fracs[o][i]
                    yerr = frac_errs[o][i]
                    plt.errorbar(x, y, yerr)

               if (
                         (row == plot_dims[0]-1 and col==plot_dims[1]-1) or
                         (i==obsids.size-1)
               ):
                    plt.tight_layout()
                    if args.pdf: pdf.savefig(fig)
                    else: plt.show()
                    plt.clf()

     if True:

          plot_dims = (5,3)
          nplots = reduce(lambda x, y: x*y, [fracs[o].shape[1] for o in fracs])

          for o in fracs:

               for i in range(fracs[o].shape[1]):

                    row = int(i/plot_dims[1]) % plot_dims[0]
                    col = i % plot_dims[1]
                    plt.subplot2grid(plot_dims, (row,col))

                    x = dates
                    y = fracs[o][:,i]
                    yerr = frac_errs[o][:,i]
                    plt.errorbar(x, y, yerr)

                    if (
                              (row == plot_dims[0]-1 and col==plot_dims[1]-1) or
                              (i==nplots-1)
                    ):
                         plt.tight_layout()
                         if args.pdf: pdf.savefig(fig)
                         else: plt.show()
                         plt.clf()

     if args.pdf:
          pdf.close()
        
     print(fracs)

if __name__ == '__main__':
    main()
