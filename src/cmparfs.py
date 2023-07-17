import astropy.io.fits
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib
matplotlib.rc('text', usetex=True)

def read_specresp(infile, args):
    hdulist = astropy.io.fits.open(infile)
    data = hdulist['specresp'].data
    energ_lo, energ_hi, specresp = data.field('energ_lo'), data.field('energ_hi'), data.field('specresp')
    hdulist.close()
    if args.energy:
        return 0.5 * (energ_lo + energ_hi), specresp
    else:
        return (0.5 * 12.39854 * (1/energ_lo + 1/energ_hi))[::-1], specresp[::-1]

def main():
    parser = argparse.ArgumentParser(
        description='Compare ARF files, plotting second/first.'
    )
    parser.add_argument('-y', '--ylabel', help='Y axis label', default=r'\textrm{Second / First}')
    parser.add_argument('-t', '--title', help='Plot title', default=r'\textrm{SPECRESP Ratio}')
    parser.add_argument('-o', '--outfile', help='Output file name')
    parser.add_argument('-e', '--energy', action='store_true', help='Plot energy rather than wavelength.')
    parser.add_argument('first')
    parser.add_argument('second')
    args = parser.parse_args()

    w1, ea1 = read_specresp(args.first, args)
    w2, ea2 = read_specresp(args.second, args)

    ratio = ea2 / np.interp(w2, w1, ea1)
    plt.plot(w2, ratio)
    if args.energy:
        plt.xlabel(r'\textrm{Energy (keV)}')
    else:
        plt.xlabel(r'\textrm{Wavelength (\AA)}')
    plt.ylabel(args.ylabel)
    plt.title(args.title)
    plt.ylim(0.9, 1.1)

    plt.tight_layout()

    if (args.outfile): plt.savefig(args.outfile)
    else: plt.show()

    plt.close()

if __name__ == '__main__':
    main()
