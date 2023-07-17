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

import hz43
import util

def main():

    parser = argparse.ArgumentParser(
        description='Determine extent of 0th order CHIP coordinates for HRC-S/LETG HZ 43 observations',
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    args = parser.parse_args()

    chipx_min, chipx_max, chipy_min, chipy_max = (None for i in range(4))

    obsids, _ = hz43.obsids_years('HRC-S')
    ci=np.linspace(0, 1, obsids.size)

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = (11, 8.5))

    plt.title(r'\textrm{HZ 43 HRC-S/LETG Zeroth Order CHIP[XY] Extents}')
    plt.xlabel(r'\textrm{CHIPY}')
    plt.ylabel(r'\textrm{CHIPX}')

    for j in range(obsids.size):
        obsid = obsids[j]

        evt2file = util.evt2_file(obsid)
        hdr = util.read_header(evt2file)

        x0, y0, r0, x, y, chipx, chipy = util.read_evt2(evt2file)

        i = np.where((x-x0)**2 + (y-y0)**2 < r0)

        if chipx_min is None:
            chipx_min = chipx[i].min()
            chipx_max = chipx[i].max()
            chipy_min = chipy[i].min()
            chipy_max = chipy[i].max()

        xmin, xmax = chipx[i].min(), chipx[i].max()
        ymin, ymax = chipy[i].min(), chipy[i].max()

        color = plt.cm.RdYlBu(ci[j])
        plt.plot([ymin,ymax], [0.5*(xmin+xmax) for i in range(2)], '-', color=color)
        plt.plot([0.5*(ymin+ymax) for i in range(2)], [xmin, xmax], '-', color=color)

        sys.stderr.write("Obsid={:d}, CHIPX=[{:d}, {:d}], CHIPY=[{:d}, {:d}]\n".format(obsid, xmin, xmax, ymin, ymax))

        if xmin < chipx_min: chipx_min = xmin
        if xmax > chipx_max: chipx_max = xmax

        if ymin < chipy_min: chipy_min = ymin
        if ymax > chipy_max: chipy_max = ymax


    sys.stderr.write("Inclusive CHIPX=[{:d}, {:d}], CHIPY=[{:d}, {:d}]\n".format(chipx_min, chipx_max, chipy_min, chipy_max))
    print("{:d} {:d} {:d} {:d}".format(chipx_min, chipx_max, chipy_min, chipy_max))

    if args.pdf:
        pdf.savefig(fig)
        pdf.close()
    else:
        plt.show()

    plt.close()

if __name__ == '__main__':
    main()
