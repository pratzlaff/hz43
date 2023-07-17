import util
import hz43
import astropy.io.fits
import numpy as np
import glob
import argparse

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc, rcParams
rc('text', usetex=True)
#rcParams.update({'font.size': 14})

def get_ratios(obsid, wav_ranges):

    d_f, h_f = util.read_pha2(util.pha2_file(obsid, tg_reprocess='tg_reprocess_pifilter'))
    d_u, h_u = util.read_pha2(util.pha2_file(obsid, tg_reprocess='tg_reprocess'))

    ratios = np.zeros((wav_ranges.shape[1],))
    ratio_errs = ratios.copy()

    for i in range(ratios.size):

        wav_min = wav_ranges[0][i]
        wav_max = wav_ranges[1][i]

        tg_m = 1
        if wav_min < 0:
            tg_m = -1

        j = np.where(d_f['tg_m']==tg_m)[0][0]

        # rate, err
        r_f, re_f = util.calc_rate(d_f['counts'][j],
                                   (d_f['background_up'][j] +
                                    d_f['background_down'][j]),
                                   h_f
                                   )
        r_u, re_u = util.calc_rate(d_u['counts'][j],
                                   (d_u['background_up'][j] +
                                    d_u['background_down'][j]),
                                   h_u
                                   )
        bin_lo = d_f['bin_lo'][j].copy()
        bin_hi = d_f['bin_hi'][j].copy()

        if wav_min < 0:
            wav_min, wav_max = -wav_max, -wav_min

        ind = np.where(np.logical_and(bin_lo>=wav_min, bin_hi<=wav_max))

        rate_f = r_f[ind].sum()
        rate_u = r_u[ind].sum()

        err_f = np.sqrt((re_f[ind]**2).sum())
        err_u = np.sqrt((re_u[ind]**2).sum())

        ratios[i] = rate_f / rate_u
        ratio_errs[i] = np.sqrt(ratios[i]**2 * ( (err_f/rate_f)**2 + (err_u/rate_u)**2))

    return ratios, ratio_errs

    
def main():

    parser = argparse.ArgumentParser(
        description='Plot filtered/unfiltered rates for HRC-S/LETG HZ 43 observations.'
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('--png', help='Output PNG file name prefix.')
    args = parser.parse_args()

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = (8.5, 11))

    obsids, dates = hz43.obsids_years('HRC-S')

    wav_ranges = np.zeros((2, 20))
    wav_ranges[0][:10] = -(np.arange(10)*10+70)[::-1]
    wav_ranges[0][10:] = (np.arange(10)*10+70)
    wav_ranges[1] = wav_ranges[0]+10

    binsize=2
    neg_min=55
    neg_max=163
    pos_min=67
    pos_max=175

    # binsize=10
    # neg_min=60
    # neg_max=160
    # pos_min=70
    # pos_max=170

    # FIXME: assumes binsize evenly divides the range
    neg_n = int((neg_max-neg_min)/binsize)
    pos_n = int((pos_max-pos_min)/binsize)

    wav_neg = -((np.arange(neg_n)+1)*binsize+neg_min)[::-1]
    wav_pos = np.arange(pos_n)*binsize+pos_min

    wav_ranges = np.zeros((2, neg_n+pos_n))
    wav_ranges[0][:neg_n] = wav_neg
    wav_ranges[0][neg_n:] = wav_pos
    wav_ranges[1] = wav_ranges[0]+binsize

    ratios = np.zeros((obsids.size, wav_ranges.shape[1]))
    ratio_errs = ratios.copy()

    for i in range(obsids.size):
        r, rerr = get_ratios(obsids[i], wav_ranges)
        ratios[i] = r
        ratio_errs[i] = rerr

    plot_dims=(3,2)
    for i in range(ratios.shape[1]):
        row = int(i/plot_dims[1]) % plot_dims[0]
        col = i % plot_dims[1]
        plt.subplot2grid(plot_dims, (row,col))

        x = dates
        y = ratios[:,i]
        err = ratio_errs[:,i]

        plt.plot(x, y, 'r.')
        #plt.errorbar(x, y, err, fmt='r.')
        plt.ylim(0.92, 1.02)

        title = r'$' + '{}'.format(wav_ranges[0][i]) + r'< \lambda < '+'{}'.format(wav_ranges[1][i]) + r'$'
        plt.title(title)

        if col == 0:
            if args.png:
                plt.ylabel(r'\textrm{Ratio}')
            else:
                plt.ylabel(r'\textrm{Rate: Filtered / Unfiltered}')
        if row == plot_dims[0]-1:
            plt.xlabel(r'\textrm{Date}')

        if (
                (row == plot_dims[0]-1 and col==plot_dims[1]-1) or
                (i==ratios.shape[1]-1)
        ):
            plt.tight_layout()
            if args.pdf:
                pdf.savefig(fig)
            elif args.png:
                 plt.savefig(args.png+'{:02d}.png'.format(1+int(i/(plot_dims[0]*plot_dims[1]))), bbox_inches='tight')
            else:
                plt.show()
            plt.clf()

    if args.pdf:
        pdf.close()
        
if __name__ == '__main__':
    main()
