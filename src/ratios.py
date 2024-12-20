import sys
import glob
import argparse
import astropy.io.fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import re
import sys

import response
import util
import hz43
import flux
import symbols

useflux=False
maxorder=None
fig=None
pdf=None
corrected=None

# unconcered about leap years

def ymd2datestr(y, m, d):
    return f'{y:04d}-{m:02d}-{d:02d}'

def ymd2frac(y, m, d):
    y = int(y)
    m = int(m)
    d = float(d)
    return y + (md2dn(m, d)-1)/365

def md2dn(m, d):
    m = int(m)
    d = float(d)
    if ( m <= 2 ):
        return (m-1)*31+d
    else:
        return int((m+1)*30.6)-63+d

def dn2md(dn):
    sum = 0
    month_days = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    for i in range(len(month_days)):
        sum += month_days[i]
        if sum >= dn:
            return i+1, dn-(sum-month_days[i])
    raise ValueError(f'{dn}: {sum}')


def frac2ymd(frac):
    y = int(frac)
    dn = int(365*(frac-y))+1
    m, d = dn2md(dn)
    return y, m, d

def datestr2ymd(datestr):
    m = re.match(r'(\d{4})-(\d{2})-(\d{2})', datestr)
    year = int(m.group(1))
    month = int(m.group(2))
    day = float(m.group(3))
    return year, month, day

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

    predicted = predicted[ind].sum()

    ratio = rate / predicted
    ratio_err = rate_err / predicted

    global useflux
    if useflux:
        ratio = flux_ / predicted
        ratio_err = flux_err / predicted

    return rate, rate_err, flux_, flux_err, ratio, ratio_err

# get HRC-S/LETG counts light curves for dispersed orders
def dispersed_lc(tg_reprocess, exclude, merge):
    global maxorder

    orders = { 'neg':-1, 'pos':+1 }

    obsids, years = hz43.obsids_years('HRC-S', exclude=exclude)
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

    if merge is not None:
        data = { obsids[i] : {'year':years[i],
                              'date_str':date_str[i],
                              'rate': { o:rates[o][:,i] for o in orders },
                              'rate_err': { o:rate_errs[o][:,i] for o in orders },
                              'flux': { o:fluxes[o][:,i] for o in orders },
                              'flux_err': { o:flux_errs[o][:,i] for o in orders },
                              'ratio': { o:ratios[o][:,i] for o in orders },
                              'ratio_err': { o:ratio_errs[o][:,i] for o in orders },
                              } for i in range(len(obsids))
                }

        merge_disp_rates(data, merge)

        obsids = np.array(list(data.keys()))
        years = np.array([ data[o]['year'] for o in obsids ])
        date_str = [ data[o]['date_str'] for o in obsids ]
        rates = { order:np.stack([data[o]['rate'][order] for o in obsids]).transpose() for order in orders }
        rate_errs = { order:np.stack([data[o]['rate_err'][order] for o in obsids]).transpose() for order in orders }
        fluxes = { order:np.stack([data[o]['flux'][order] for o in obsids]).transpose() for order in orders }
        flux_errs = { order:np.stack([data[o]['flux_err'][order] for o in obsids]).transpose() for order in orders }
        ratios = { order:np.stack([data[o]['ratio'][order] for o in obsids]).transpose() for order in orders }
        ratio_errs = { order:np.stack([data[o]['ratio_err'][order] for o in obsids]).transpose() for order in orders }

    return obsids, years, date_str, w1, w2, rates, rate_errs, fluxes, flux_errs, ratios, ratio_errs

def merge_disp_rates(data, merge):
    orders = ('neg', 'pos')
    for a in merge:
        try:

            rates = { order:np.array([ data[o]['rate'][order] for o in a ]) for order in orders }
            rate_errs = { order:np.array([ data[o]['rate_err'][order] for o in a ]) for order in orders }
            w = { order:1/rate_errs[order]/rate_errs[order] for order in orders }
            rate = { order:(rates[order] * w[order]).sum(axis=0) / w[order].sum(axis=0) for order in orders }
            rate_err = { order:np.sqrt(1/w[order].sum(axis=0)) for order in orders }

            fluxes = { order:np.array([ data[o]['flux'][order] for o in a ]) for order in orders }
            flux_errs = { order:np.array([ data[o]['flux_err'][order] for o in a ]) for order in orders }
            w = { order:1/flux_errs[order]/flux_errs[order] for order in orders }
            flux = { order:(fluxes[order] * w[order]).sum(axis=0) / w[order].sum(axis=0) for order in orders }
            flux_err = { order:np.sqrt(1/w[order].sum(axis=0)) for order in orders }

            # FIXME: is this valid for the ratios?
            ratios = { order:np.array([ data[o]['ratio'][order] for o in a ]) for order in orders }
            ratio_errs = { order:np.array([ data[o]['ratio_err'][order] for o in a ]) for order in orders }
            w = { order:1/ratio_errs[order]/ratio_errs[order] for order in orders }
            ratio = { order:(ratios[order] * w[order]).sum(axis=0) / w[order].sum(axis=0) for order in orders }
            ratio_err = { order:np.sqrt(1/w[order].sum(axis=0)) for order in orders }

            data[a[-1]]['rate'] = rate
            data[a[-1]]['rate_err'] = rate_err

            data[a[-1]]['flux'] = flux
            data[a[-1]]['flux_err'] = flux_err

            data[a[-1]]['ratio'] = ratio
            data[a[-1]]['ratio_err'] = ratio_err

            for i in range(len(a)-1):
                del data[a[i]]
        except:
            raise

# get HRC-S/LETG counts light curves for zeroth order
def zeroth_lc(detector, tg_reprocess, exclude, merge):
    if (detector == 'HRC-S'):
        obsids, years = hz43.obsids_years('HRC-S', exclude=exclude)
    elif (detector == 'HRC-I'):
        obsids, years = hz43.obsids_years('HRC-I', exclude=exclude)
    else:
        raise ValueError(det)

    rates, rate_errs, exposures = util.zeroth_rates(obsids, tg_reprocess=tg_reprocess)
    model_rates = hz43.predicted_rates(obsids)

    if merge is not None:
        data = { obsids[i] : {'year':years[i],
                              'rate':rates[i],
                              'rate_err':rate_errs[i],
                              'exposure':exposures[i],
                              'model_rate':model_rates[i],
                              } for i in range(len(obsids))
                }

        merge_zero_rates(data, merge)

        obsids = np.array(list(data.keys()))
        years = np.array([ data[o]['year'] for o in data ])
        rates = np.array([ data[o]['rate'] for o in data ])
        rate_errs = np.array([ data[o]['rate_err'] for o in data ])
        model_rates = np.array([ data[o]['model_rate'] for o in data ])
        exposures = np.array([ data[o]['exposure'] for o in data ])

    date_obs = []
    for i in range(len(obsids)):
        date_obs.append(util.read_header(util.pha2_file(obsids[i], tg_reprocess=tg_reprocess))['date-obs'][0:10])

    return obsids, years, date_obs, rates, rate_errs, model_rates, rates/model_rates, rate_errs/model_rates

def merge_zero_rates(data, merge):

    # each element is an array of obsids to merge
    for a in merge:
        try:
            rates = np.array([data[o]['rate'] for o in a])
            rate_errs = np.array([data[o]['rate_err'] for o in a])
            rate_errs = np.array([data[o]['rate_err'] for o in a])
            model_rates = np.array([data[o]['model_rate'] for o in a])
            exposures = np.array([data[o]['exposure'] for o in a])
            w = 1 /rate_errs / rate_errs
            rate = (rates * w).sum() / w.sum()
            rate_err = np.sqrt(1/w.sum())

            model_rate = (model_rates*exposures).sum() / exposures.sum()

            data[a[-1]]['rate'] = rate
            data[a[-1]]['rate_err'] = rate_err
            data[a[-1]]['model_rate'] = model_rate

            for i in range(len(a)-1):
                del data[a[i]]
        except:
            pass

def plot_zero(d, args, label=None, relative=True):
    x = d['year']
    if args.absolute:
        y = d['rate']
        yerr = d['rate_err']
    else:
        y = d['ratio']
        yerr = d['ratio_err']
    #print(x, y, yerr)
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

def plot_disp_wavdep_ratios(lc_disp, args):
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
            ratio_err = lc_disp['ratio_err'][order][1:][:,i]
            plt.errorbar(wav, ratio, ratio_err, fmt=fmt[order])
            plt.title(f"{lc_disp['obsid'][i]} - {lc_disp['date-obs'][i]}")

        if args.ymin is not None:
            plt.ylim(bottom=args.ymin)
        if args.ymax is not None:
            plt.ylim(top=args.ymax)

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

def collect_data(lc_0, lc_disp):
    d = {
        'obsid' : lc_0['obsid'],
        'year' : lc_0['year'],
        'date-obs' : lc_0['date-obs'],
        'r_0' : lc_0['ratio'],
        'rerr_0' : lc_0['ratio_err'],
        'r_pos' : lc_disp['ratio']['pos'][0],
        'rerr_pos' : lc_disp['ratio_err']['pos'][0],
        'r_neg' : lc_disp['ratio']['neg'][0],
        'rerr_neg' : lc_disp['ratio_err']['neg'][0],
    }
    return d

# get HRC-S QEU file specifications
def qeu_params(lc_0, lc_disp, args):

    # return values will be
    # ---------------------
    # cvsd - self-explanatory
    # r0 - ratio of observed to predicted for zeroth order
    # rpos - ratio of observed to predicted for positive orders
    # rneg - ratio of observed to predicted for negative orders
    # obsid - obsid whose residual ratios are to be used for wavelength-dependent corrections

    cvsd = ['1999-07-22']
    for i in range(10):
        cvsd.append(f'{2000+i:04d}-01-01')
    obsid = [None] * len(cvsd)
    cvsd_year = [ymd2frac(*datestr2ymd(d)) for d in cvsd]

    d = collect_data(lc_0, lc_disp)
    ind = np.where(d['year'] > 2010)[0]

    hv_1_date = ymd2frac(2012, 3, 29)
    hv_2_date = ymd2frac(2021, 5, 14)
    hv_3_date = ymd2frac(2024, 9, 20)
    hv_changes = { 14422 : '2012-03-29',
                   24575 : '2021-05-14',
                   78427 : '2024-09-20',
                 }
    for i in ind:

        o = d['obsid'][i]

        # FIXME: special case
        if o == 78427:
            continue
        if o == 28427:
            cvsd_year.append(0.5*(d['year'][i]+d['year'][i-2]))
        else:
            cvsd_year.append(0.5*(d['year'][i]+d['year'][i-1]))
        cvsd.append(ymd2datestr(*frac2ymd(cvsd_year[-1])))

        obsid.append(o)
        if o in hv_changes:
            cvsd[-1] = hv_changes[o]
            cvsd_year[-1] = ymd2frac(*datestr2ymd(cvsd[-1]))

        # FIXME: special case
        if o == 29452:
            obsid.append(78427)
            cvsd.append(hv_changes[78427])
            cvsd_year.append(ymd2frac(*datestr2ymd(cvsd[-1])))

    # these are the dates of the middle of the time period each file
    # will cover
    cvsd_year = np.array(cvsd_year)
    year_eff = np.copy(cvsd_year)
    year_eff[:-1] = 0.5 * (cvsd_year[1:] + cvsd_year[:-1])

    r0 = np.zeros(len(cvsd))
    rpos = r0.copy()
    rneg = r0.copy()

    # fit ratios before first HV change
    ind = np.where(d['year'] < hv_1_date)[0]
    p_0 = np.polynomial.Polynomial.fit(d['year'][ind],
                                        d['r_0'][ind],
                                        1,
                                        w=1/d['rerr_0'][ind]
                                        )
    p_pos = np.polynomial.Polynomial.fit(d['year'][ind],
                                          d['r_pos'][ind],
                                          1,
                                          w=1/d['rerr_pos'][ind]
                                        )
    p_neg = np.polynomial.Polynomial.fit(d['year'][ind],
                                          d['r_neg'][ind],
                                          1,
                                          w=1/d['rerr_neg'][ind]
                                        )
    sys.stderr.write('p_0: '+np.array2string(p_0.convert().coef) + '\n')
    sys.stderr.write('p_pos: '+np.array2string(p_pos.convert().coef) + '\n')
    sys.stderr.write('p_neg: '+np.array2string(p_neg.convert().coef) + '\n')
    ind = np.where(year_eff < hv_1_date)[0]
    r0[ind] = p_0(year_eff[ind])
    rpos[ind] = p_pos(year_eff[ind])
    rneg[ind] = p_neg(year_eff[ind])

    # fit ratios between first and second HV changes
    ind = np.where((d['year'] > hv_1_date) & (d['year'] < hv_2_date))[0]
    p_0 = np.polynomial.Polynomial.fit(d['year'][ind],
                                        d['r_0'][ind],
                                        1,
                                        w=1/d['rerr_0'][ind]
                                        )
    p_pos = np.polynomial.Polynomial.fit(d['year'][ind],
                                          d['r_pos'][ind],
                                          1,
                                          w=1/d['rerr_pos'][ind]
                                        )
    p_neg = np.polynomial.Polynomial.fit(d['year'][ind],
                                          d['r_neg'][ind],
                                          1,
                                          w=1/d['rerr_neg'][ind]
                                        )
    sys.stderr.write('p_0: '+np.array2string(p_0.convert().coef) + '\n')
    sys.stderr.write('p_pos: '+np.array2string(p_pos.convert().coef) + '\n')
    sys.stderr.write('p_neg: '+np.array2string(p_neg.convert().coef) + '\n')
    ind = np.where((year_eff > hv_1_date) & (year_eff < hv_2_date))[0]
    r0[ind] = p_0(year_eff[ind])
    rpos[ind] = p_pos(year_eff[ind])
    rneg[ind] = p_neg(year_eff[ind])

    # just duplicate observed ratios after 2nd HV change
    ind = np.where(year_eff > hv_2_date)[0]
    in_n = d['obsid'].size
    out_n = len(cvsd)

    for i in ind:
        j = -(out_n-i)

        # FIXME: special cases resulting from 78427
        if i == 35:
            j -= 2
        if i == 33 or i == 34:
            j += 1

        if obsid[i] != d['obsid'][j]:
            raise ValueError(f'{i}\t{j}\t{in_n}\t{out_n}\t{obsid[i]}\t{d["obsid"][j]}')
        r0[i] = d['r_0'][j]
        rpos[i] = d['r_pos'][j]
        rneg[i] = d['r_neg'][j]

    print('\t'.join(('obsid', 'cvsd', 'year_eff', 'r0', 'rpos', 'rneg')))
    print('\t'.join(('S',)*2+('N',)*4))
    for i in range(len(cvsd)):
        print('\t'.join((str(obsid[i]), cvsd[i], f'{year_eff[i]:.5f}', f'{r0[i]:.3f}', f'{rpos[i]:.3f}', f'{rneg[i]:.3f}')))

    return cvsd, obsid, r0, rpos, rneg

def write_disp_ratios(lc_0, lc_disp):
    d = collect_data(lc_0, lc_disp)
    for i in range(len(lc_0['obsid'])):
        year, month, day = frac2ymd(d['year'][i])
        ymd = f'{year:04d}-{month:02d}-{day:02d}'
        sys.stderr.write('\t'.join((
            str(d['obsid'][i]),
            d['date-obs'][i],
            str(d['year'][i]),
            ymd,
            f'{d["r_0"][i]:.3f}',
            f'{d["rerr_0"][i]:.4f}',
            f'{d["r_pos"][i]:.3f}',
            f'{d["rerr_pos"][i]:.4f}',
            f'{d["r_neg"][i]:.3f}',
            f'{d["rerr_neg"][i]:.4f}',
        )) + '\n')

def lc_0(detnam, tg_reprocess, exclude, merge):
    lc_0 = {}
    lc_0.update(zip(('obsid', 'year', 'date-obs', 'rate', 'rate_err', 'model_rate', 'ratio', 'ratio_err'), zeroth_lc(detnam, tg_reprocess, exclude, merge)))
    return lc_0

def lc_disp(tg_reprocess, exclude, merge):
    lc_disp = {}
    lc_disp.update(zip(('obsid', 'year', 'date-obs', 'bin_lo', 'bin_hi', 'rate', 'rate_err', 'flux', 'flux_err', 'ratio', 'ratio_err'), dispersed_lc(tg_reprocess, exclude, merge)))
    return lc_disp


def select_file_index(cvsd, cvsds):
    for i in range(len(cvsds)):
        if cvsds[i] > cvsd:
            return i-1
    return i

def qeu_correct(lc_0, lc_disp, args):
    # special cases of HV changes
    hv0 = (14324, 14396, 14397) # taken 2012-07-04
    hv1 = ( 14238, )             # taken 2012-03-18

    cvsd, obsid, r0, rpos, rneg = qeu_params(lc_0, lc_disp, args)
    for i in range(len(lc_0['obsid'])):
        date_obs = lc_0['date-obs'][i]
        obsid = lc_0['obsid'][i]

        if obsid in hv0:
            date_obs='2012-03-16'
        if obsid in hv1:
            date_obs='2012-03-30'
        if obsid==62635:
            date_obs='2021-05-15'

        j = select_file_index(lc_0['date-obs'][i], cvsd)

        lc_0['ratio'][i] /= r0[j]
        lc_0['ratio_err'][i] /= r0[j]

        lc_disp['ratio']['pos'][:,i] /= rpos[j]
        lc_disp['ratio_err']['pos'][:,i] /= rpos[j]

        lc_disp['ratio']['neg'][:,i] /= rneg[j]
        lc_disp['ratio_err']['neg'][:,i] /= rneg[j]

def mkplots(i_0, s_0, disp, args):
    global fig, pdf, corrected, useflux

    if useflux:
        ylabel = 'Flux / Predicted'
    else:
        ylabel = 'Rate / Predicted'
    if not args.noi:
        plot_zero(i_0, args, label='HRC-I: 0th')
    if not args.nos:
        plot_zero(s_0, args, label=r'HRC-S: 0th')
        for order in disp['rate']:
            plot_dispersed(disp, order, 0)
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
    parser.add_argument('-a', '--absolute', help='Absolute rates.', action='store_true')
    parser.add_argument('-m', '--maxorder', help='Maximum ARF/RMF order to read.', default=3, type=int)
    #parser.add_argument('-e','--exclude', nargs='*', type=int, default=[24958,62635,25615], help='Exclude obsids')
    parser.add_argument('-e','--exclude', nargs='*', type=int, help='Exclude obsids')
    parser.add_argument('--ymin', type=float, help='Lower Y plot limit.')
    parser.add_argument('--ymax', type=float, help='Upper Y plot limit.')
    parser.add_argument('--merge', action='append', nargs='+', type=int)

    args = parser.parse_args()

    global maxorder, useflux

    useflux = args.flux
    maxorder = args.maxorder

    hrci_lc_0 = None
    hrcs_lc_0 = None
    hrcs_lc_disp = None

    if not args.noi:
        hrci_lc_0 = lc_0('HRC-I', args.tg_reprocess_hrci, args.exclude, args.merge)

    if not args.nos:
        hrcs_lc_0 = lc_0('HRC-S', args.tg_reprocess_hrcs, args.exclude, args.merge)
        hrcs_lc_disp = lc_disp(args.tg_reprocess_hrcs, args.exclude, args.merge)

    figsize = (11, 8.5)

    global fig, pdf
    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = figsize)

    mkplots(hrci_lc_0, hrcs_lc_0, hrcs_lc_disp, args)

    if args.nos:
        if args.pdf:
            pdf.close()
        sys.exit()

    plot_disp_wavdep_ratios(hrcs_lc_disp, args)
    write_disp_ratios(hrcs_lc_0, hrcs_lc_disp)

    if args.corrected:
        qeu_correct(hrcs_lc_0, hrcs_lc_disp, args)
        mkplots(hrci_lc_0, hrcs_lc_0, hrcs_lc_disp, args)
        plot_disp_wavdep_ratios(hrcs_lc_disp, args)
        write_disp_ratios(hrcs_lc_0, hrcs_lc_disp)

    if args.resdir:
        write_residuals(hrcs_lc_disp, args.resdir)

    if args.pdf:
        pdf.close()

    sys.exit()

if __name__ == '__main__':
    main()
