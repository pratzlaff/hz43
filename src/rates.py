import sys
import glob
import argparse
import astropy.io.fits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import hz43
import util
import symbols

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

# get HRC-S/LETG rates for dispersed orders
def dispersed_rates(tg_reprocess='tg_reprocess', merge=None):
    orders = { 'neg':-1, 'pos':+1 }

    obsids, years = hz43.obsids_years('HRC-S')
    w1, w2 = wav_ranges()

    date_str = []

    rates, rate_errs = ({} for i in range(2))
    for order in w1:
        rates[order] = np.zeros((w1[order].size, years.size))
        rate_errs[order] = rates[order].copy()

    for i in range(obsids.size):
        obsid = obsids[i]

        # read PHA2
        pha2  = util.pha2_file(obsid, tg_reprocess=tg_reprocess)
        d, h = util.read_pha2(util.pha2_file(obsid, tg_reprocess=tg_reprocess))
        date_str.append(h['date-obs'][:10])

        rows = {}
        for order in w1:
            ind = np.where(d['tg_m']==orders[order])[0][0]
            src = d['counts'][ind]
            bg = d['background_up'][ind] + d['background_down'][ind]
            bin_lo = d['bin_lo'][ind]
            bin_hi = d['bin_hi'][ind]
            
            for j in range(w1[order].size):
                ind2 = np.where((bin_lo>=w1[order][j]) & (bin_hi<w2[order][j]))
                rates[order][j][i], rate_errs[order][j][i] = util.calc_rate(src[ind2].sum(), bg[ind2].sum(), h)

    data = { obsids[i] : {'year':years[i],
                          'date_str':date_str[i],
                          'rate': { 'neg':rates['neg'][:,i],
                                    'pos':rates['pos'][:,i],
                                   },
                          'rate_err': { 'neg':rate_errs['neg'][:,i],
                                        'pos':rate_errs['pos'][:,i],
                                   },
                          } for i in range(len(obsids))
            }

    merge_disp_rates(data, merge)

    obsids = list(data.keys())
    years = np.array([ data[o]['year'] for o in obsids ])
    date_str = [ data[o]['date_str'] for o in obsids ]
    rates = { 'neg':np.stack([data[o]['rate']['neg'] for o in obsids]).transpose(),
              'pos':np.stack([data[o]['rate']['pos'] for o in obsids]).transpose(),
             }
    rates = { order:np.stack([data[o]['rate'][order] for o in obsids]).transpose() for order in orders }
    rate_errs = { order:np.stack([data[o]['rate_err'][order] for o in obsids]).transpose() for order in orders }

    return obsids, years, date_str, w1, w2, rates, rate_errs

def merge_disp_rates(data, merge):
    # each element is an array of obsids to merge
    orders = ('neg', 'pos')
    for a in merge:
        try:
            rates = { order:np.array([ data[o]['rate'][order] for o in a ]) for order in orders }
            rate_errs = { order:np.array([ data[o]['rate_err'][order] for o in a ]) for order in orders }
            w = { order:1/rate_errs[order]/rate_errs[order] for order in orders }
            rate = { order:(rates[order] * w[order]).sum(axis=0) / w[order].sum(axis=0) for order in orders }
            rate_err = { order:np.sqrt(1/w[order].sum(axis=0)) for order in orders }
            data[a[-1]]['rate'] = rate
            data[a[-1]]['rate_err'] = rate_err
            for i in range(len(a)-1):
                del data[a[i]]
        except:
            raise

# get HRC rates for 0th order
def zeroth_rates(detector, tg_reprocess='tg_reprocess', merge=None):
    if (detector == 'HRC-S'):
        obsids, years = hz43.obsids_years('HRC-S')
    elif (detector == 'HRC-I'):
        obsids, years = hz43.obsids_years('HRC-I')
    else:
        raise ValueError(det)

    rates, rate_errs, exposures = util.zeroth_rates(obsids, tg_reprocess=tg_reprocess)

    data = { obsids[i] : {'year':years[i],
                          'rate':rates[i],
                          'rate_err':rate_errs[i],
                          } for i in range(len(obsids))
            }

    merge_zero_rates(data, merge)

    years = np.array([ data[o]['year'] for o in data ])
    rates = np.array([ data[o]['rate'] for o in data ])
    rate_errs = np.array([ data[o]['rate_err'] for o in data ])

    return years, rates, rate_errs

def merge_zero_rates(data, merge):

    # each element is an array of obsids to merge
    for a in merge:
        try:
            rates = np.array([data[o]['rate'] for o in a])
            rate_errs = np.array([data[o]['rate_err'] for o in a])
            w = 1 /rate_errs / rate_errs
            rate = (rates * w).sum() / w.sum()
            rate_err = np.sqrt(1/w.sum())
            data[a[-1]]['rate'] = rate
            data[a[-1]]['rate_err'] = rate_err
            for i in range(len(a)-1):
                del data[a[i]]
        except:
            pass

def plot_0th(args, detector):
    tg_reprocess = { 'HRC-I' : args.tg_reprocess_hrci,
                     'HRC-S' : args.tg_reprocess_hrcs,
                    }.get(detector)
    rates_0 = {}
    rates_0.update(zip(('year', 'rate', 'rate_err'), zeroth_rates(detector, tg_reprocess=tg_reprocess, merge=args.merge)))
    year = rates_0['year']
    rate = rates_0['rate']
    rate_err = rates_0['rate_err']
    ratio = rate / rate[0]
    ratio_err = np.sqrt( ratio**2 * ( (rate_err/rate)**2 + (rate_err[0]/rate[0])**2 ) )
    if args.absolute:
        y = rate
        yerr = rate_err
        ylabel = 'Rate'
    else:
        y = ratio
        yerr = ratio_err
        ylabel = 'Rate Ratio'
    plt.errorbar(year, y, yerr, label=f'{detector}: 0th')
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel(ylabel)

def plot_disp(args, rates_disp):
    d=rates_disp
    labels = { 'pos' : 'HRC-S: +1st', 'neg' : 'HRC-S: -1st' }
    linestyles = { 'pos' : '-', 'neg' : '--' }
    year = d['year']
    for order in d['bin_lo']:
        rate = d['rate'][order][0,:]
        rate_cmp = rate[0]
        rate_err = d['rate_err'][order][0,:]
        rate_err_cmp = rate_err[0]
        ratio = rate/rate_cmp
        ratio_err = ratio * np.sqrt( (rate_err/rate)**2 + (rate_err_cmp/rate_cmp)**2 )

        y = ratio
        yerr = ratio_err
        if args.absolute:
            y = rate
            yerr = rate_err
        
        label = labels[order] + f": {d['bin_lo'][order][0]:.0f}-{d['bin_hi'][order][0]:.0f} {symbols.ANGSTROM}"
        plt.errorbar(year, y, yerr, label=label, linestyle=linestyles[order])

def plot_disp_wavdep(args, rates_disp):
    global fig, pdf
    fmt = { 'pos' : 'k-', 'neg' : 'r-' }
    plot_dims = (2, 3)
    for i in range(rates_disp['year'].size):
        row = int(i/plot_dims[1]) % plot_dims[0]
        col = i % plot_dims[1]

        plt.subplot2grid(plot_dims, (row,col))
        for order in rates_disp['bin_lo']:
            wav = 0.5*(rates_disp['bin_lo'][order][1:] + rates_disp['bin_hi'][order][1:])
            rate = rates_disp['rate'][order][1:][:,i]
            rate_cmp = rates_disp['rate'][order][1:][:,0]
            rate_err = rates_disp['rate_err'][order][1:][:,i]
            rate_err_cmp = rates_disp['rate_err'][order][1:][:,0]
            ratio = rate/rate_cmp
            ratio_err = ratio * np.sqrt( (rate_err/rate)**2 + (rate_err_cmp/rate_cmp)**2 )

            y = ratio
            yerr = ratio_err
            ylabel = 'Rate Ratio'
            if args.absolute:
                y = rate
                yerr = rate_err
                ylabel = 'Rate'
                
            plt.errorbar(wav, y, yerr, fmt=fmt[order])
            plt.title(f"{rates_disp['obsid'][i]} - {rates_disp['date'][i]}")

        plt.grid()
        if args.ymin is not None:
            plt.ylim(bottom=args.ymin)
        if args.ymax is not None:
            plt.ylim(top=args.ymax)

        if (row==plot_dims[0]-1) or (i>=rates_disp['year'].size-plot_dims[1]):
            plt.xlabel(f'{symbols.LAMBDA} ({symbols.ANGSTROM})')
        if col==0:
            plt.ylabel(ylabel)

        if (
            (row == plot_dims[0]-1 and col==plot_dims[1]-1) or
                (i==rates_disp['year'].size-1)
        ):
            plt.tight_layout()
            if args.pdf:
                pdf.savefig(fig)
            else:
                plt.show()
            plt.clf()

def main():
    global fig, pdf
    parser = argparse.ArgumentParser(
        description='Plot HRC HZ 43 rate ratios, compared to the first HZ 43 observation',
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('--tg_reprocess_hrci', default='tg_reprocess', help='tg_reprocess output directory for HRC-I.')
    parser.add_argument('--tg_reprocess_hrcs', default='tg_reprocess', help='tg_reprocess output directory for HRC-S.')
    parser.add_argument('--noi', help='Do not plot I curves.', action='store_true')
    parser.add_argument('--nos', help='Do not plot S curves.', action='store_true')
    parser.add_argument('--nodisp', help='Do not plot individual observation S dispersed rates.', action='store_true')
    parser.add_argument('-a', '--absolute', help='Plot rates rather than ratios.', action='store_true')
    parser.add_argument('--ymin', type=float, help='Lower Y plot limit.')
    parser.add_argument('--ymax', type=float, help='Upper Y plot limit.')
    parser.add_argument('-m', '--merge', type=int, action='append', nargs='+', default=[[25615,27916],[25614,29077]])
    parser.add_argument('--width', type=float, default=11, help='PDF width in inches.')
    parser.add_argument('--height', type=float, default=8.5, help='PDF height in inches.')
    parser.add_argument('--lw', type=float, default=1, help='Line widths.')
    parser.add_argument('--fs', type=float, default=10, help='Font sizes.')
    args = parser.parse_args()

    matplotlib.rcParams['lines.linewidth'] = args.lw
    matplotlib.rcParams['axes.linewidth'] = args.lw
    matplotlib.rcParams['font.size'] = args.fs

    if args.pdf:
        pdf = PdfPages(args.pdf)
        figsize = (args.width, args.height)
        plot_dims = (1, 1)
        fig = plt.figure(figsize=figsize)

    if not args.noi:
        plot_0th(args, 'HRC-I')

    if not args.nos:
        plot_0th(args, 'HRC-S')
        s_disp = {}
        s_disp.update(zip(('obsid', 'year', 'date', 'bin_lo', 'bin_hi', 'rate', 'rate_err'), dispersed_rates(tg_reprocess=args.tg_reprocess_hrcs, merge=args.merge)))
        plot_disp(args, s_disp)


    title = 'HRC HZ 43 Count Rate Ratios'
    if args.absolute:
        title = 'HRC HZ 43 Count Rates'
    plt.title(title)
    plt.legend()
    plt.tight_layout()

    if args.pdf:
        pdf.savefig(fig)
    else:
        plt.show()
    plt.clf()

    if not args.nos and not args.nodisp:
        plot_disp_wavdep(args, s_disp)

    if args.pdf:
        pdf.close()

    sys.exit()

if __name__ == '__main__':
    main()
