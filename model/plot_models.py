import glob
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib
matplotlib.rc('text', usetex=True)

def plot_model(year, order, xlabel=True):

    ostr = { -1:'neg', 1:'pos' }[order]

    model_file = year + '_' + ostr + '.model'
    wavm, ratiom = np.loadtxt(model_file, unpack=True)

    obs_file = 'ratios/' + year + '_' + ostr + '.ratio'
    wav, flux, ratio, err = np.loadtxt(obs_file, unpack=True)

    plt.plot(wavm, ratiom, 'k-')
    plt.plot(wav, ratio, 'k.')
    plt.errorbar(wav, ratio, err, ecolor='k', fmt='none')

    plt.ylim(0.85, 1.1)

    plt.title(r'$\textrm{%s - TG\_M=%s}$' % ( year, { -1:'-1', 1:'+1'}[order] ))

    if xlabel: plt.xlabel(r'$\textrm{Wavelength (\AA)}$')

    ylabel = { -1:r'$\textrm{Flux ratio}$', 1:'' }[order]
    plt.ylabel({ -1:r'$\textrm{Flux ratio}$', 1:'' }[order])

def main():
    parser = argparse.ArgumentParser(
        description='Plot CARMA models of HZ 43 flux ratios'
    )
    args = parser.parse_args()

    files = glob.glob('2[0-9][0-9][0-9]_*.model')
    years = { }

    for f in files:
        years[f[:4]] = None

    years = years.keys()
    years.sort()

    plot_dims = (3,2)

    for i in xrange(len(years)):
        year = years[i]
        row = i % plot_dims[0]

        end_page = not (i+1) % plot_dims[0] or i==len(years)-1

        plt.subplot2grid(plot_dims, (row, 0))
        plot_model(year, -1, xlabel=end_page)
        plt.subplot2grid(plot_dims, (row, 1))
        plot_model(year, +1, xlabel=end_page)

        if end_page:
            plt.tight_layout()
            plt.savefig('model_plots_%d.png' % int(i/plot_dims[0]), bbox_inches='tight')
            plt.show()



if __name__ == '__main__':
    main()
