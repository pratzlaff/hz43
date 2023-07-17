import sys
import glob
import argparse
import astropy.io.fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)
#rcParams.update({'font.size': 14})

import response
import util
import hz43
import flux

# These are the dates just prior to a new piecemeal linear fit of flux
# vs time
ichanges=(2011.14, 2015.13, 2018.28, 2021.1284)
schanges=(2012.3, 2021.2217)

useflux=False
maxorder=None

def wav_ranges():

    # wavelength intervals will be (e.g.):
    # TG_M = -1 : 60-160, 60-65, 65-70, ..., 155-160
    # TG_M = +1 : 70-170, 70-75, 75-80, ..., 165-170

    # FIXME: starting with N0011, obsid 1011 negative-order residuals
    # are not being calculated, and wavelength range doesn't seem to
    # be the culprit

    increment = { 'neg' : 2, 'pos' : 2 }
    steps = { 'neg' : 50, 'pos' : 52 }
    minimum = { 'neg' : 57, 'pos' : 69 }

    w1 = { order : (np.arange(steps[order]+1)-1)*increment[order]+minimum[order]
           for order in increment }
    sys.stderr.write(str(w1)+"\n")

    w2 = {}
    for order in w1:
        w2[order] = w1[order] + increment[order]
        w1[order][0] = w1[order][1]
        w2[order][0] = w2[order][-1]

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

def plot_zero(d, ratio=False, label=None, relative=True):
    x = d['year']

    if ratio:
        y = d['ratio']
        yerr = d['ratio_err']
    else:
        y = d['rate']
        yerr = d['rate_err']
        if relative:
            y_2008_5 = np.interp(2008.5, x, y)
            y = y/y_2008_5
            yerr = yerr/y_2008_5

    plt.errorbar(x, y, yerr, label=label)

def plot_dispersed(d, order, index, ratio=False, color=None):
    global useflux
    labels = { 'pos' : 'HRC-S: +1st', 'neg' : 'HRC-S: -1st' }

    linestyles = { 'pos' : '-', 'neg' : '--' }

    x = d['year']

    if ratio is False:
        y = d['rate'][order][index]
        yerr = d['rate_err'][order][index]
        y_2008_5 = np.interp(2008.5, x, y)
        y = y/y_2008_5
        yerr = yerr/y_2008_5
    else:
        y = d['ratio'][order][index]
        yerr = d['ratio_err'][order][index]

    label = labels[order] + ": {:.0f}-{:.0f} \\AA".format(d['bin_lo'][order][index], d['bin_hi'][order][index])
    plt.errorbar(x, y, yerr, label=label, color=color, linestyle=linestyles[order])

def x_1(p):
    return (1-p[1])/p[0]

def fit_s(x, y_, yerr_):

    y_2008_5 = np.interp(2008.5, x, y_)
    y_2011 = np.interp(2010.8, x, y_)

    # y = y_/y_2008_5
    # yerr = yerr_/y_2008_5

    y=y_
    yerr=yerr_

    i1 = np.where(x<2012)[0]
    p1, cov1 = np.polyfit(x[i1], y[i1], 1, w=yerr[i1], cov=True)

    # y = y_/y_2011
    # yerr = yerr_/y_2011

    i2 = np.where((x>2012) & (x<2021.2217) & (x==x) & (y==y) & (yerr==yerr))[0]
    sys.stderr.write(str(x[i2]) + str(y[i2]) + str(yerr[i2]) + "\n")
    p2, cov2 = np.polyfit(x[i2], y[i2], 1, w=yerr[i2], cov=True)

    i3 = np.where((x>2021.2217) & (x==x) & (y==y) & (yerr==yerr))[0]
    sys.stderr.write(str(x[i3]) + str(y[i3]) + str(yerr[i3]) + "\n")
    p3, cov3 = np.polyfit(x[i3], y[i3], 1, w=yerr[i3], cov=True)

    sigma1 = np.sqrt(np.diagonal(cov1))
    sigma2 = np.sqrt(np.diagonal(cov2))
    sigma3 = np.sqrt(np.diagonal(cov3))

    return p1, sigma1, p2, sigma2, p3, sigma3

def fit_i(x, y_, yerr_):

    y_2008_5 = np.interp(2008.5, x, y_)
    y_2011 = np.interp(2010.8, x, y_)

    # y = y_/y_2008_5
    # yerr = yerr_/y_2008_5

    y=y_
    yerr=yerr_

    # FIXME: need to nail down the date change
    i1 = np.where(x<2021)[0]
    p1, cov1 = np.polyfit(x[i1], y[i1], 1, w=yerr[i1], cov=True)

    i2 = np.where(x>2021)[0]
    p2, cov2 = np.polyfit(x[i2], y[i2], 1, w=yerr[i2], cov=True)

    sigma1 = np.sqrt(np.diagonal(cov1))
    sigma2 = np.sqrt(np.diagonal(cov2))

    return p1, sigma1, p2, sigma2

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

def plot_higher_order_contributions(args):
    
    return

    global maxorder

    model_wav, model_flux = hz43.model()

    obsids, _ = hz43.obsids_years('HRC-S')

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = (11, 8.5))

    plt.title(r'\textrm{HZ 43 HRC-S/LETG}')
    plt.xlabel(r'\textrm{Date}')
    plt.xlabel(r'$\lambda$ \textrm{(\AA)}')

    bin_lo, bin_hi, rmfs = util.get_rmfs('HRC-S', maxorder)
    wav = 0.5*(bin_lo+bin_hi)
    energy = util.w2e(wav)

    # for j in range(obsids.size):
    for j in range(1):
        obsid = obsids[j]
        pha2 = util.pha2_file(obsid, tg_reproces=tg_reprocess_hrcs)
        hdr = util.read_header(pha2)

        jnk, jnk, garfs = util.get_garfs(obsid, maxorder)
        resp = { order : garfs[order] * rmfs[order] for order in garfs }
        rate_predicted = { order : np.interp(wav, model_wav, model_flux) * (bin_hi-bin_lo) / 1.602e-9 / energy * resp[order] for order in resp }

        sys.stderr.write(resp['neg'].shape + "\n")

        maxi = rate_predicted['pos'].shape[0]
        if 0:
            plt.ylabel(r'\textrm{log rate}')
            plt.yscale('log')
            for i in range(maxi):
                plt.plot(wav[0], rate_predicted['pos'][i], label='{:+d}'.format(i+1))
                plt.plot(wav[0], rate_predicted['neg'][i], label='{:+d}'.format(-i-1))

        elif 1:
            plt.ylabel(r'\textrm{fractional order contributions')
            for i in range(1,maxi):
                plt.plot(wav[0], rate_predicted['pos'][i] / rate_predicted['pos'][0:maxi].sum(axis=0), label='{:+d}'.format(i+1))
                plt.plot(wav[0], rate_predicted['neg'][i] / rate_predicted['neg'][0:maxi].sum(axis=0), label='{:+d}'.format(-i-1))

    plt.legend()

    if args.pdf:
        pdf.savefig(fig)
        pdf.close()
    else:
        plt.show()
        
    sys.exit()

def main():
    parser = argparse.ArgumentParser(
        description='HRC-S/LETG HZ 43 observations',
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('-r', '--resdir', help='Save residual ratios to named directory.')
    parser.add_argument('--tg_reprocess_hrci', default='tg_reprocess', help='tg_reprocess output directory for HRC-I.')
    parser.add_argument('--tg_reprocess_hrcs', default='tg_reprocess', help='tg_reprocess output directory for HRC-S.')
    parser.add_argument('-f', '--flux', help='Compare with model fluxes, rather than rates.', action='store_true')
    parser.add_argument('-c', '--corrected', help='Plot corrected ratio curves.', action='store_true')
    parser.add_argument('--rate', help='Plot rate curves and exit. May be combined with --pred.', action='store_true')
    parser.add_argument('--pred', help='Plot predicted flux curves and exit. May be combined with --rate.', action='store_true')
    parser.add_argument('--noi', help='Do not plot I curves.', action='store_true')
    parser.add_argument('--nos', help='Do not plot S curves.', action='store_true')
    parser.add_argument('-m', '--maxorder', help='Maximum ARF/RMF order to read.', default=3, type=int)

    args = parser.parse_args()

    global maxorder, useflux

    useflux = args.flux
    maxorder = args.maxorder

    plot_higher_order_contributions(args)

    if not args.noi:
        hrci_lc_0 = {}
        hrci_lc_0.update(zip(('year', 'rate', 'rate_err', 'model_rate', 'ratio', 'ratio_err'), zeroth_lc('HRC-I', tg_reprocess=args.tg_reprocess_hrci)))

    if not args.nos:
        hrcs_lc_0 = {}
        hrcs_lc_0.update(zip(('year', 'rate', 'rate_err', 'model_rate', 'ratio', 'ratio_err'), zeroth_lc('HRC-S', tg_reprocess=args.tg_reprocess_hrcs)))

        lc_1 = {}
        lc_1.update(zip(('obsid', 'year', 'date', 'bin_lo', 'bin_hi', 'rate', 'rate_err', 'flux', 'flux_err', 'ratio', 'ratio_err'), dispersed_lc(tg_reprocess=args.tg_reprocess_hrcs)))

    if False:
        plot_zero(hrci_lc_0, label=r'HRC-S: 0th', relative=False)
        plt.ylabel('Rate')
        plt.xlabel('Date')
        plt.title('HZ 43, HRC-I/LETG, 0th Order')
        plt.savefig('hz43_hrci_0th_rate.png')
        plt.clf()
        exit(0)

    plot_dims = (3, 1)
    figsize = (8.5, 11)

    if args.rate + args.pred == 1:
        rcParams.update({'font.size': 14})
        figsize = (11, 8.5)
        plot_dims = (1, 1)
    elif args.rate + args.pred == 2:
        plot_dims = (2, 1)

    if args.pdf:
        rcParams.update({'font.size': 14})
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = figsize)


    if args.rate:
        plt.subplot2grid(plot_dims, (0,0))
        if not args.noi:
            plot_zero(hrci_lc_0, label='HRC-I: 0th')
        if not args.nos:
            plot_zero(hrcs_lc_0, label=r'HRC-S: 0th')
            print(hrcs_lc_0['rate'])
            for order in lc_1['rate']:
                plot_dispersed(lc_1, order, 0)
                print(lc_1['rate'][order][0])
        plt.legend(fontsize=10)
        plt.title(r'\textrm{HZ 43: HRC/LETG Count Rates}')
        plt.ylabel(r'\textrm{Rate / 2008.5}')
        if not args.pred:
            plt.xlabel(r'\textrm{Year}')

    if args.pred:
        if args.rate:
            plt.subplot2grid(plot_dims, (1,0))
        else:
            plt.subplot2grid(plot_dims, (0,0))
        if useflux:
            ylabel = r'\textrm{Flux / Predicted}'
        else:
            ylabel = r'\textrm{Rate / Predicted}'
        if not args.noi:
            plot_zero(hrci_lc_0, label='HRC-I: 0th', ratio=True)
        if not args.nos:
            plot_zero(hrcs_lc_0, label=r'HRC-S: 0th', ratio=True)
            for order in lc_1['rate']:
                plot_dispersed(lc_1, order, 0, ratio=True)
        if not args.rate:
            plt.legend(fontsize=10)
            plt.title(r'\textrm{HZ 43: HRC/LETG Ratios to Predicted}')
        plt.ylabel(ylabel)
        plt.xlabel(r'\textrm{Year}')

    if args.rate or args.pred:
        plt.tight_layout()
        if args.pdf:
            pdf.savefig(fig)
            pdf.close()
        else:
            plt.show()
        exit(0)

    #
    # plot rate light curves for 0th order and outer plates
    #

    ylabel = r'\textrm{Rate / 2008.5}'
    plt.subplot2grid(plot_dims, (0,0))
    if not args.noi:
        plot_zero(hrci_lc_0, label='HRC-I: 0th')
    if not args.nos:
        plot_zero(hrcs_lc_0, label=r'HRC-S: 0th')
        for order in lc_1['rate']:
            plot_dispersed(lc_1, order, 0)
    plt.legend(fontsize=10)
    plt.title(r'\textrm{HZ 43: HRC/LETG Count Rates}')
    plt.ylabel(ylabel)

    #
    # plot rates light curves for each wavelength bin separately
    #

    if not args.nos:
        ci = { order : np.linspace(1, 0, lc_1['bin_lo'][order].size-1) for order in lc_1['bin_lo'] }

        rows = { 'pos' : 1, 'neg' : 2 }
        titles = { 'pos' : r'\textrm{Positive Orders}', 'neg' : r'\textrm{Negative Orders}' } 

        for order in lc_1['bin_lo']:
            plt.subplot2grid(plot_dims, (rows[order],0))
            for i in range(1, lc_1['bin_lo'][order].size):
                plot_dispersed(lc_1, order, i, color=plt.cm.RdYlBu(ci[order][i-1]))
            plt.xlabel(r'\textrm{Date}')
            plt.ylabel(ylabel)
            plt.title(titles[order])

    plt.tight_layout()

    if args.pdf:
        pdf.savefig(fig)
    else:
        plt.show()

    #
    # plot ratio light curves, [rate|flux] / model, for 0th order and outer plates
    #

    if useflux:
        ylabel = r'\textrm{Flux / Predicted}'
    else:
        ylabel = r'\textrm{Rate / Predicted}'
    plt.subplot2grid(plot_dims, (0,0))
    if not args.noi:
        plot_zero(hrci_lc_0, label='HRC-I: 0th', ratio=True)
    if not args.nos:
        plot_zero(hrcs_lc_0, label=r'HRC-S: 0th', ratio=True)
        for order in lc_1['rate']:
            plot_dispersed(lc_1, order, 0, ratio=True)
    plt.title(r'\textrm{HZ 43: HRC/LETG Ratios to Predicted}')
    plt.ylabel(ylabel)

    if args.nos:
        plt.tight_layout()

        if args.pdf:
            pdf.savefig(fig)
        else:
            plt.show()

        sys.exit()

    #
    # plot ratio light curves, [rate|flux] / model, for each wavelength bin separately
    #

    for order in lc_1['bin_lo']:
        plt.subplot2grid(plot_dims, (rows[order],0))
        for i in range(1, lc_1['bin_lo'][order].size):
            plot_dispersed(lc_1, order, i, color=plt.cm.RdYlBu(ci[order][i-1]), ratio=True)
        plt.xlabel(r'\textrm{Date}')
        plt.ylabel(ylabel)
        plt.title(titles[order])

    plt.tight_layout()

    if args.pdf:
        pdf.savefig(fig)
    else:
        plt.show()

    if args.corrected:

        # fitting
        if not args.noi:
            i_p1_0th, i_err1_0th, i_p2_0th, i_err2_0th = fit_i(hrci_lc_0['year'], hrci_lc_0['ratio'], hrci_lc_0['ratio_err'])
            sys.stderr.write('HRC-I 0th' + "\n")
            sys.stderr.write('\tHV0: slope={:.3g} ({:.3g}), xoff={:.6g}'.format(i_p1_0th[0], i_err1_0th[0], x_1(i_p1_0th)) + "\n")
            sys.stderr.write('\tHV1: slope={:.3g} ({:.3g}), xoff={:.6g}'.format(i_p2_0th[0], i_err2_0th[0], x_1(i_p2_0th)) + "\n")
            sys.stderr.write('' + "\n")

        s_p1_0th, s_err1_0th, \
        s_p2_0th, s_err2_0th, \
        s_p3_0th, s_err3_0th \
        = fit_s(hrcs_lc_0['year'], hrcs_lc_0['ratio'], hrcs_lc_0['ratio_err'])

        s_p1_neg, s_err1_neg, s_p2_neg, s_err2_neg, s_p3_neg, s_err3_neg = fit_s(lc_1['year'], lc_1['ratio']['neg'][0], lc_1['ratio_err']['neg'][0])
        s_p1_pos, s_err1_pos, s_p2_pos, s_err2_pos, s_p3_pos, s_err3_pos = fit_s(lc_1['year'], lc_1['ratio']['pos'][0], lc_1['ratio_err']['pos'][0])
        sys.stderr.write("HRC-S 0th\n")
        sys.stderr.write(f"\tHV0: "
                         f"slope={s_p1_0th[0]:.3g} ({s_err1_0th[0]:.3g}), "
                         f"xoff={x_1(s_p1_0th):.6g}\n"
                         )
        sys.stderr.write(f"\tHV1: "
                         f"slope={s_p2_0th[0]:.3g} ({s_err2_0th[0]:.3g}), "
                         f"xoff={x_1(s_p2_0th):.6g}\n"
                         )
        sys.stderr.write(f"\tHV2: "
                         f"slope={s_p3_0th[0]:.3g} ({s_err3_0th[0]:.3g}), "
                         f"xoff={x_1(s_p3_0th):.6g}\n"
                         )
        sys.stderr.write("\n")

        sys.stderr.write("HRC-S -1\n")
        sys.stderr.write(f"\tHV0: "
                         f"slope={s_p1_neg[0]:.3g} ({s_err1_neg[0]:.3g}), "
                         f"xoff={x_1(s_p1_neg):.6g}\n"
                         )
        sys.stderr.write(f"\tHV1: "
                         f"slope={s_p2_neg[0]:.3g} ({s_err2_neg[0]:.3g}), "
                         f"xoff={x_1(s_p2_neg):.6g}\n"
                         )
        # sys.stderr.write(f"\tHV2: "
        #                  f"slope={s_p3_neg[0]:.3g} ({s_err3_neg[0]:.3g}), "
        #                  f"xoff={x_1(s_p3_neg):.6g}\n"
        #                  )
        sys.stderr.write("\n")

        sys.stderr.write("HRC-S +1\n")
        sys.stderr.write(f"\tHV0: "
                         f"slope={s_p1_pos[0]:.3g} ({s_err1_pos[0]:.3g}), "
                         f"xoff={x_1(s_p1_pos):.6g}\n"
                         )
        sys.stderr.write(f"\tHV1: "
                         f"slope={s_p2_pos[0]:.3g} ({s_err2_pos[0]:.3g}), "
                         f"xoff={x_1(s_p2_pos):.6g}\n"
                         )
        # sys.stderr.write(f"\tHV2: "
        #                  f"slope={s_p3_pos[0]:.3g} ({s_err3_pos[0]:.3g}), "
        #                  f"xoff={x_1(s_p3_pos):.6g}\n"
        #                  )
        sys.stderr.write("\n")

        s_slope_hv0_0th = s_p1_0th[0]
        s_slope_hv1_0th = s_p2_0th[0]
        #s_slope_hv2_0th = s_p3_0th[0]

        s_slope_hv0_0th_err = s_err1_0th[0]
        s_slope_hv1_0th_err = s_err2_0th[0]
        #s_slope_hv2_0th_err = s_err3_0th[0]

        s_slope_hv0_1st = 0.5 * (s_p1_neg[0]+s_p1_pos[0])
        s_slope_hv1_1st = 0.5 * (s_p2_neg[0]+s_p2_pos[0])
        #s_slope_hv2_1st = 0.5 * (s_p3_neg[0]+s_p3_pos[0])

        s_slope_hv0_1st_err = 0.5 * np.sqrt(s_err1_neg[0]**2+s_err1_pos[0]**2)
        s_slope_hv1_1st_err = 0.5 * np.sqrt(s_err2_neg[0]**2+s_err2_pos[0]**2)
        #s_slope_hv2_1st_err = 0.5 * np.sqrt(s_err3_neg[0]**2+s_err3_pos[0]**2)

        s_xoff_hv0 = 2008.5
        s_xoff_hv1 = 1/3. * (x_1(s_p2_0th) + x_1(s_p2_neg) + x_1(s_p2_pos)) # 2010.9
        #s_xoff_hv2 = 1/3. * (x_1(s_p3_0th) + x_1(s_p3_neg) + x_1(s_p3_pos))

        sys.stderr.write("HRC-S HV_0: "
                         f"xoff={s_xoff_hv0:.6g}, "
                         f"0th slope={s_slope_hv0_0th:.3g} ({s_slope_hv0_0th_err:.3g}), "
                         f"1st slope={s_slope_hv0_1st:.3g} ({s_slope_hv0_1st_err:.3g})\n"
                         )

        sys.stderr.write("HRC-S HV_1: "
                         f"xoff={s_xoff_hv1:.6g}, "
                         f"0th slope={s_slope_hv1_0th:.3g} ({s_slope_hv1_0th_err:.3g}), "
                         f"1st slope={s_slope_hv1_1st:.3g} ({s_slope_hv1_1st_err:.3g})\n"
                         )
        # sys.stderr.write("HRC-S HV_2: "
        #                  f"xoff={s_xoff_hv2:.6g}, "
        #                  f"0th slope={s_slope_hv2_0th:.3g} ({s_slope_hv2_0th_err:.3g}), "
        #                  f"1st slope={s_slope_hv2_1st:.3g} ({s_slope_hv2_1st_err:.3g})\n"
        #                  )
        sys.stderr.write("\n")

        print("{:.6g}\t{:.3g}\t{:.3g}\t{:.6g}\t{:.3g}\t{:.3g}".format(
            s_xoff_hv0, s_slope_hv0_0th, s_slope_hv0_1st,
            s_xoff_hv1, s_slope_hv1_0th, s_slope_hv1_1st,
            #s_xoff_hv2, s_slope_hv2_0th, s_slope_hv2_1st
        ))

        #
        # 0th order corrections
        #
        
        correction = np.ones(hrcs_lc_0['year'].size)

        index = np.where(hrcs_lc_0['year']<2012)[0]
        correction[index] = (s_xoff_hv0 - hrcs_lc_0['year'][index]) * s_slope_hv0_0th

        # the 2nd HV change was on 2021-05-14 (day 134), but for HZ43
        # obsid discrimination purposes, the dividing line can be
        # placed at 2021.2217, which was the date separating the previous
        # HV observation of HZ 43, and the test of the new HV
        year_hv2 = 2021+133/365.
        year_hv2 = 2021.2217

        index = np.where((hrcs_lc_0['year']>2012) & (hrcs_lc_0['year']<year_hv2))[0]
        correction[index] = (s_xoff_hv1 - hrcs_lc_0['year'][index]) * s_slope_hv1_0th

        # grey decline frozen in place at the 2nd HV change
        index = np.where(hrcs_lc_0['year'] > year_hv2)[0]
        correction[index] = (s_xoff_hv1 - year_hv2) * s_slope_hv1_0th

        print('0th order corrections: ', correction, '\n')
        hrcs_lc_0['ratio'] /= (1-correction)
        hrcs_lc_0['ratio_err'] /= (1-correction)

        #
        # 1st order corrections
        #

        correction = np.zeros(lc_1['year'].size)

        index = np.where(lc_1['year']<2012)[0]
        correction[index] = (s_xoff_hv0 - lc_1['year'][index]) * s_slope_hv0_1st

        index = np.where((lc_1['year']>2012) & (lc_1['year']<year_hv2))[0]
        correction[index] = (s_xoff_hv1 - lc_1['year'][index]) * s_slope_hv1_1st

        index = np.where(lc_1['year']>year_hv2)[0]
        correction[index] = (s_xoff_hv1 - year_hv2) * s_slope_hv1_1st

        print('dispersed order corrections: ', correction, '\n')
        for order in lc_1['ratio']:
            lc_1['ratio'][order] /= (1-correction)
            lc_1['ratio_err'][order] /= (1-correction)

        #
        # plot ratio light curves, [rate|flux] / model, for 0th order and outer plates
        #

        if useflux:
            ylabel = r'\textrm{Flux / Predicted}'
        else:
            ylabel = r'\textrm{Rate / Predicted}'

        plt.subplot2grid(plot_dims, (0,0))
        if not args.noi:
            plot_zero(hrci_lc_0, label='HRC-I: 0th', ratio=True)

        if not args.nos:
            plot_zero(hrcs_lc_0, label=r'HRC-S: 0th', ratio=True)
            for order in lc_1['rate']:
                plot_dispersed(lc_1, order, 0, ratio=True)
            plt.title(r'\textrm{HZ 43: HRC/LETG Ratios to Predicted}')
            plt.ylabel(ylabel)

            #
            # plot ratio light curves, [rate|flux] / model, for each wavelength bin separately
            #
            for order in lc_1['bin_lo']:
                plt.subplot2grid(plot_dims, (rows[order],0))
                for i in range(1, lc_1['bin_lo'][order].size):
                    plot_dispersed(lc_1, order, i, color=plt.cm.RdYlBu(ci[order][i-1]), ratio=True)
            plt.xlabel(r'\textrm{Date}')
            plt.ylabel(ylabel)

        plt.tight_layout()

        if args.pdf:
            pdf.savefig(fig)
        else:
            plt.show()

        # end if args.corrected

    #
    # ratio plots
    #
    
    plot_dims = (3, 2)

    for i in range(lc_1['year'].size):

        row = int(i/plot_dims[1]) % plot_dims[0]
        col = i % plot_dims[1]

        plt.subplot2grid(plot_dims, (row,col))
        plot_dispersed_ratio(lc_1, i)

        if row == plot_dims[0]-1:
            plt.xlabel(r'$\lambda$ \textrm{(\AA)}')
        if col==0:
            if useflux:
                ylabel  = r'\textrm{Flux / Predicted}'
            else:
                ylabel  = r'\textrm{Rate / Predicted}'
            plt.ylabel(ylabel)

        if (
                (row == plot_dims[0]-1 and col==plot_dims[1]-1) or
                (i==lc_1['year'].size-1)
        ):
            plt.tight_layout()
            if args.pdf:
                pdf.savefig(fig)
            else:
                plt.show()
            plt.clf()

    if args.pdf:
        pdf.close()

    if args.resdir:
        write_residuals(lc_1, args.resdir)

    sys.exit()

if __name__ == '__main__':
    wav_ranges()
    main()
