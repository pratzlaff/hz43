import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import sys

import response

def main():
    parser = argparse.ArgumentParser(
        description='Print ARF values on a wavelength grid.'
    )
    parser.add_argument('arffile', help='ARF filename.')

    args = parser.parse_args()
    wlo, whi, ea = response.read_arf(args.arffile)
    wav = 0.5*(wlo+whi)
    plt.plot(wav, ea)
    plt.show()

if __name__ == '__main__':
    main()
