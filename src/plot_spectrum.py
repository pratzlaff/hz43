import argparse
import util
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)
#rcParams.update({'font.size': 14})

np.seterr(divide='ignore', invalid='ignore')

def main():
    parser = argparse.ArgumentParser(
        description='Plot a spectrum.'
    )
    parser.add_argument('-c', '--combine', type=int, default=80, help='Number of channels to combine for each bin.')
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('obsid', help='Obsid to Plot.')
    args = parser.parse_args()

    spectrum = util.get_spectrum(args.obsid, 'LEG', args.combine)

    if args.pdf:
        pdf = PdfPages(args.pdf)

    fig = plt.figure(figsize = (11, 8.5))

    labels = { 'neg' : r'\textrm{TG\_M=-1}', 'pos' : r'\textrm{TG\_M=+1}' }
    symbols = { 'neg' : 'r--', 'pos' : 'b-' }
    for order in spectrum:

        plt.errorbar(spectrum[order]['lambda'],
                     spectrum[order]['rate'],
                     spectrum[order]['rate_err'],
                     fmt=symbols[order],
                     label=labels[order])

    plt.xlabel(r'$\lambda$ (\AA)')
    plt.xlim(60,180)
    plt.ylim(1e-2)
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.tight_layout()

    if args.pdf:
        pdf.savefig(fig)
        pdf.close()
    else:
        plt.show()

    plt.close()

if __name__ == '__main__':
    main()
