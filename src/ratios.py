import sys
import glob
import argparse
import astropy.io.fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import response
import util
import hz43
import flux
import symbols

useflux=False
maxorder=None
fig=None
pdf=None

def wav_ranges():

    # wavelength intervals will be (e.g.):
    # TG_M = -1 : 57-157, 57-59, 59-61, ..., 155-157
    # TG_M = +1 : 69-173, 69-71, 71-73, ..., 171-173

    increment = { 'neg' : 2, 'pos' : 2 }
    steps = { 'neg' : 50, 'pos' : 52 }
    minimum = { 'neg' : 57, 'pos' : 69 }

    w1 = { order : (np.arange(steps[order]+1)-1)*increment[order]+minimum[order]
           for order in increment }
    w2 = {}
    for order in w1:
        w2[order] = w1[order] + increment[order]
        w1[order][0] = w1[order][1]
        w2[order][0] = w2[order][-1]

    sys.stderr.write(str(w1)+"\n")
    sys.stderr.write(str(w2)+"\n")

    return w1, w2


def dispersed_flux(src, bg, resp, bin_lo, bin_hi, wav_lo, wav_hi, hdr, predicted, factor):
    ind = np.where((bin_lo>=wav_lo) & (bin_hi<wav_hi))
    rate, rate_err, flux_, flux_err = flux.flux_summed(src[ind], bg[ind], resp[ind], bin_lo[ind], bin_hi[ind], hdr, factor[ind])

    ratio = rate / predicted[ind].sum()
    ratio_err = rate_err / predicted[ind].sum()

    global useflux
    if useflux:
        ratio = flux_ / predicted[ind].sum()
        ratio_err = flux_err / predicted[ind].sum()

    return rate, rate_err, flux_, flux_err, ratio, ratio_err

# get HRC-S/LETG counts light curves for dispersed orders
def dispersed_lc(tg_reprocess='tg_reprocess'):
    global maxorder

    orders = { 'neg':-1, 'pos':+1 }

    obsids, years = hz43.obsids_years('HRC-S')
    w1, w2 = wav_ranges()

    model_flux = None

    date_str = []

    rates, rate_errs, fluxes, flux_errs, ratios, ratio_errs = ({} for  i in range(6))
    for order in w1:
        rates[order] = np.zeros((w1[order].size, years.size))
        rate_errs[order] = rates[order].copy()
        fluxes[order] = rates[order].copy()
        flux_errs[order] = rates[order].copy()
        ratios[order] = rates[order].copy()
        ratio_errs[order] = rates[order].copy()

    for i in range(obsids.size):
        obsid = obsids[i]

        bin_lo, bin_hi, resp = response.get_response(obsid, 'LEG', maxorder=maxorder)

        if model_flux is None:
            model_wav, model_flux = hz43.model()
            wav = 0.5*(bin_lo+bin_hi)
            energy = util.w2e(wav)
            model_flux = np.interp(wav, model_wav, model_flux) * (bin_hi - bin_lo)

        # ratio of 1st-order counts / all orders
        predicted_rate = { order : model_flux / 1.602e-9 / energy * resp[order] for order in resp }
        maxi = predicted_rate['pos'].shape[0]
        factor = { order : predicted_rate[order][0] / predicted_rate[order].sum(axis=0) for order in resp }

        # read PHA2
        d, h = util.read_pha2(util.pha2_file(obsid, tg_reprocess=tg_reprocess))

        date_str.append(h['date-obs'][:10])

        rows = {}
        for order in w1:

            ind = np.where(d['tg_m']==orders[order])[0][0]
            src = d['counts'][ind]
            bg = d['background_up'][ind] + d['background_down'][ind]
            
            predicted = predicted_rate[order][0]
            global useflux
            if useflux:
                predicted = model_flux[0]
            
            for j in range(w1[order].size):
                r, rerr, f, ferr, ra, raerr = dispersed_flux(src, bg, resp[order][0], bin_lo[0], bin_hi[0], w1[order][j], w2[order][j], h, predicted, factor[order])
                
                rates[order][j][i] = r
                rate_errs[order][j][i] = rerr
                fluxes[order][j][i] = f
                flux_errs[order][j][i] = ferr
                ratios[order][j][i] = ra
                ratio_errs[order][j][i] = raerr

    return obsids, years, date_str, w1, w2, rates, rate_errs, fluxes, flux_errs, ratios, ratio_errs

# get HRC-S/LETG counts light curves for zeroth order
def zeroth_lc(detector, tg_reprocess='tg_reprocess'):
    if (detector == 'HRC-S'):
        obsids, years = hz43.obsids_years('HRC-S')
    elif (detector == 'HRC-I'):
        obsids, years = hz43.obsids_years('HRC-I')
    else:
        raise ValueError(det)

    rates, rate_errs = util.zeroth_rates(obsids, tg_reprocess=tg_reprocess)
    model_rates = hz43.predicted_rates(obsids)
    return years, rates, rate_errs, model_rates, rates/model_rates, rate_errs/model_rates

def plot_zero(d, label=None, relative=True):
    x = d['year']
    y = d['ratio']
    yerr = d['ratio_err']
    plt.errorbar(x, y, yerr, label=label)

def plot_dispersed(d, order, index, color=None):
    global useflux
    labels = { 'pos' : 'HRC-S: +1st', 'neg' : 'HRC-S: -1st' }
    linestyles = { 'pos' : '-', 'neg' : '--' }
    x = d['year']
    y = d['ratio'][order][index]
    yerr = d['ratio_err'][order][index]
    label = labels[order] + f": {d['bin_lo'][order][index]:.0f}-{d['bin_hi'][order][index]:.0f} {symbols.ANGSTROM}"
    plt.errorbar(x, y, yerr, label=label, color=color, linestyle=linestyles[order])

def plot_dispersed_ratio(lc, i):

    order = 'neg'
    x = 0.5*(lc['bin_lo'][order][1:] + lc['bin_hi'][order][1:])
    y = lc['ratio'][order][1:,i]
    yerr = lc['ratio_err'][order][1:,i]
    plt.errorbar(x, y, yerr, fmt='r-')

    order = 'pos'
    x = 0.5*(lc['bin_lo'][order][1:] + lc['bin_hi'][order][1:])
    y = lc['ratio'][order][1:,i]
    yerr = lc['ratio_err'][order][1:,i]
    plt.errorbar(x, y, yerr, fmt='k-')

    plt.title('{} - {}'.format(lc['obsid'][i], lc['date'][i]))

def write_residuals(lc, resdir):
    util.mkdir_p(resdir)
    for i in range(lc['obsid'].size):
        obsid = lc['obsid'][i]
        for order in lc['ratio']:
            filename = resdir + '/{}_{}.ratio'.format(obsid, order)
            wav = 0.5*(lc['bin_lo'][order][1:] + lc['bin_hi'][order][1:])
            ratio = lc['ratio'][order][1:,i]
            ratio_err = lc['ratio_err'][order][1:,i]
            np.savetxt(filename, np.transpose([wav, ratio, ratio_err]), fmt=['%.1f', '%.5f', '%.5f'], delimiter="\t", header="lambda\tflux_ratio\tflux_ratio_err")

def plot_disp_wavdep_ratios(args, lc_disp):
    global fig, pdf
    fmt = { 'pos' : 'k-', 'neg' : 'r-' }
    plot_dims = (2, 3)
    for i in range(lc_disp['year'].size):
        row = int(i/plot_dims[1]) % plot_dims[0]
        col = i % plot_dims[1]

        plt.subplot2grid(plot_dims, (row,col))
        for order in lc_disp['bin_lo']:
            wav = 0.5*(lc_disp['bin_lo'][order][1:] + lc_disp['bin_hi'][order][1:])
            ratio = lc_disp['ratio'][order][1:][:,i]
            ratio_err = lc_disp['rate_err'][order][1:][:,i]
            plt.errorbar(wav, ratio, ratio_err, fmt=fmt[order])
            plt.title(f"{lc_disp['obsid'][i]} - {lc_disp['date'][i]}")

        if (row==plot_dims[0]-1) or (i>=lc_disp['year'].size-plot_dims[1]):
            plt.xlabel(f'{symbols.LAMBDA} ({symbols.ANGSTROM})')
        if col==0:
            ylabel = 'Observed/Predicted'
            plt.ylabel(ylabel)

        if (
            (row == plot_dims[0]-1 and col==plot_dims[1]-1) or
                (i==lc_disp['year'].size-1)
        ):
            plt.tight_layout()
            if args.pdf:
                pdf.savefig(fig)
            else:
                plt.show()
            plt.clf()

def lc_0(detnam, tg_reprocess):
    lc_0 = {}
    lc_0.update(zip(('year', 'rate', 'rate_err', 'model_rate', 'ratio', 'ratio_err'), zeroth_lc(detnam, tg_reprocess)))
    return lc_0

def lc_disp(tg_reprocess):
    lc_disp = {}
    lc_disp.update(zip(('obsid', 'year', 'date', 'bin_lo', 'bin_hi', 'rate', 'rate_err', 'flux', 'flux_err', 'ratio', 'ratio_err'), dispersed_lc(tg_reprocess=tg_reprocess)))
    return lc_disp

def main():
    parser = argparse.ArgumentParser(
        description='Plot ratios of observed/predicted rates for HRC HZ 43 observations',
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('-r', '--resdir', help='Save residual ratios to named directory.')
    parser.add_argument('--tg_reprocess_hrci', default='tg_reprocess', help='tg_reprocess output directory for HRC-I.')
    parser.add_argument('--tg_reprocess_hrcs', default='tg_reprocess', help='tg_reprocess output directory for HRC-S.')
    parser.add_argument('-f', '--flux', help='Compare with model fluxes, rather than rates.', action='store_true')
    parser.add_argument('-c', '--corrected', help='Plot corrected ratio curves.', action='store_true')
    parser.add_argument('--noi', help='Do not plot I curves.', action='store_true')
    parser.add_argument('--nos', help='Do not plot S curves.', action='store_true')
    parser.add_argument('-m', '--maxorder', help='Maximum ARF/RMF order to read.', default=3, type=int)

    args = parser.parse_args()

    global maxorder, useflux

    useflux = args.flux
    maxorder = args.maxorder

    if not args.noi:
        hrci_lc_0 = lc_0('HRC-I', args.tg_reprocess_hrci)

    if not args.nos:
        hrcs_lc_0 = lc_0('HRC-S', args.tg_reprocess_hrcs)
        hrcs_lc_disp = lc_disp(args.tg_reprocess_hrcs)

    figsize = (11, 8.5)

    global fig, pdf
    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = figsize)

    #
    # plot ratio light curves, [rate|flux] / model, for 0th order and outer plates
    #
    if useflux:
        ylabel = 'Flux / Predicted'
    else:
        ylabel = 'Rate / Predicted'
    if not args.noi:
        plot_zero(hrci_lc_0, label='HRC-I: 0th')
    if not args.nos:
        plot_zero(hrcs_lc_0, label=r'HRC-S: 0th')
        for order in hrcs_lc_disp['rate']:
            plot_dispersed(hrcs_lc_disp, order, 0)
    plt.title('HZ 43: HRC/LETG Ratios to Predicted')
    plt.ylabel(ylabel)
    plt.xlabel('Year')
    plt.legend()
    plt.tight_layout()
    if args.pdf:
        pdf.savefig(fig)
    else:
        plt.show()
    plt.clf()

    if args.nos:
        if args.pdf:
            pdf.close()
        sys.exit()

    if not args.nos:
        plot_disp_wavdep_ratios(args, hrcs_lc_disp)

    if args.pdf:
        pdf.close()

    if args.resdir:
        write_residuals(hrcs_lc_disp, args.resdir)

    sys.exit()

if __name__ == '__main__':
    main()
