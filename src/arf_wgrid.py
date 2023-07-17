import argparse
import numpy as np
from scipy import interpolate
import sys

import response

def main():
    parser = argparse.ArgumentParser(
        description='Print ARF values on a wavelength grid.'
    )
    parser.add_argument('arffile', help='ARF filename.')
    parser.add_argument('wmin', type=float, help='Minimum wavelength.')
    parser.add_argument('wmax', type=float, help='Maximum wavelength.')
    parser.add_argument('n', type=int, help='Number of wavelengths.')

    args = parser.parse_args()

    wav = np.arange(args.n)*(args.wmax-args.wmin)/(args.n-1)+args.wmin

    wlo, whi, ea = response.read_arf(args.arffile)
    f = interpolate.interp1d(0.5*(wlo+whi), ea, fill_value='extrapolate')
    np.savetxt(sys.stdout, np.stack((wav, f(wav)), axis=1), fmt="%.2f\t%.6g")

if __name__ == '__main__':
    main()
