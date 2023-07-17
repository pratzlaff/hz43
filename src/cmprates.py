import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)
rcParams.update({'font.size': 14})

def read_vinay(vfile):
    obsid, rate, err = np.loadtxt(vfile, skiprows=1, usecols=(0,8,9), unpack=True);
    return obsid, rate/0.9, err/0.9

def read_pete(pfile):
    obsid, date, rate, err = np.loadtxt(pfile, usecols=(0,1,2,3), unpack=True)
    return obsid, date, rate, err

def main():

    parser = argparse.ArgumentParser(
        description='Plot ratios of HRC-I observed rates, Pete / Vinay',
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('vfile', help="Vinay's rates.")
    parser.add_argument('pfile', help="Pete's rates.")
    args = parser.parse_args()

    o_v, r_v, e_v = read_vinay(args.vfile)
    o_p, d_p, r_p, e_p = read_pete(args.pfile)

    dates = []
    ratios = []
    errs = []
    for i in range(o_p.size):
        w = np.where(o_p[i] == o_v)
        if len(w):
            k = w[0][0]
            dates.append(d_p[i])
            ratios.append(r_p[i] / r_v[k])
            errs.append(np.sqrt(ratios[-1]**2 * ( (e_p[i]/r_p[i])**2 + (e_v[k]/r_v[k])**2)))
    dates = np.array(dates)
    ratios = np.array(ratios)
    errs = np.array(errs)

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = (11, 8.5))

    plt.errorbar(dates, ratios, errs, fmt='bo')

    plt.title(r'\textrm{HZ 43: HRC-I 0th Order Observed Rates Ratios}')
    plt.xlabel(r'\textrm{Year}')
    plt.ylabel(r'\textrm{Rates: Pete / Vinay}')

    plt.tight_layout()

    if args.pdf:
        pdf.savefig(fig)
        pdf.close()
    else:
        plt.show()

    plt.close()


if __name__ == '__main__':
    main()
