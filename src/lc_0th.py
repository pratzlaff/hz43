import util
import hz43
import glob
import argparse
import astropy.io.fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)
rcParams.update({'font.size': 14})

def main():

    parser = argparse.ArgumentParser(
        description='Plot observed HRC/LETG zeroth-order count rates vs HZ 43 model',
    )
    parser.add_argument('--nohrci', action='store_true', help='Do not plot HRC-I.')
    parser.add_argument('--nohrcs', action='store_true', help='Do not plot HRC-S.')
    parser.add_argument('--printem', action='store_true', help='Print obsids, rates, errors.')
    parser.add_argument('--tg_reprocess', default='tg_reprocess', help='tg_reprocess output directory.')
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    args = parser.parse_args()

    if (args.nohrci and args.nohrcs):
        raise ValueError("Both --nohrc[is] options were given.")

    model_wav, model_flux = hz43.model()

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = (11, 8.5))

    if not args.nohrcs:
        s_obsid, s_year = hz43.obsids_years('HRC-S')
        s_predicted = hz43.predicted_rates(s_obsid)
        s_rate, s_err = util.zeroth_rates(s_obsid, tg_reprocess=args.tg_reprocess)
        s_ratio = s_rate / s_predicted
        s_ratio_err = s_err / s_predicted
        plt.errorbar(s_year, s_ratio, s_ratio_err, fmt='ro', label=r'\textrm{HRC-S}')
        if args.printem:
            for i in range(len(s_obsid)):
                print('{}\t{}\t{}\t{}'.format(s_obsid[i], s_year[i], s_rate[i], s_err[i]))

    if not args.nohrci:
        i_obsid, i_year = hz43.obsids_years('HRC-I')
        i_predicted = hz43.predicted_rates(i_obsid)
        i_rate, i_err = util.zeroth_rates(i_obsid)
        i_ratio = i_rate / i_predicted
        i_ratio_err = i_err / i_predicted
        plt.errorbar(i_year, i_ratio, i_ratio_err, fmt='bs', label=r'\textrm{HRC-I}')
        if args.printem:
            for i in range(len(i_obsid)):
                print('{}\t{}\t{}\t{}'.format(i_obsid[i], i_year[i], i_rate[i], i_err[i]))

    plt.title(r'\textrm{HZ 43: HRC/LETG Zeroth Order Light Curves}')
    plt.xlabel(r'\textrm{Year}')
    plt.ylabel(r'\textrm{Counts: Observed / Predicted}')

    plt.legend(loc='lower left')

    plt.tight_layout()

    if args.pdf:
        pdf.savefig(fig)
        pdf.close()
    else:
        plt.show()

    plt.close()


if __name__ == '__main__':
    main()
