import argparse
import numpy as np
import astropy.io.fits

def e2w(energy):
    return 12.39854/energy

def w2e(wavelength):
    return e2w(wavelength)

def hz43_model():
    return np.loadtxt('/data/legs/rpete/flight/hz43_model/non_LTE/15.model', unpack=True)

def read_arf(filename):
    hdulist = astropy.io.fits.open(filename)
    header = hdulist['specresp'].header
    data = hdulist['specresp'].data
    hdulist.close()

    energ_lo = data.field('energ_lo')[::-1]
    energ_hi = data.field('energ_hi')[::-1]
    specresp = data.field('specresp')[::-1]

    return e2w(energ_hi), e2w(energ_lo), specresp

def main():
    parser = argparse.ArgumentParser(
        description='Calculate expected HZ 43 rate from a given ARF'
    )
    parser.add_argument('infile', help='Input ARF.')
    args = parser.parse_args()

    bin_lo, bin_hi, specresp = read_arf(args.infile)
    wav = 0.5 * (bin_lo + bin_hi)
    model_wav, model_flux = hz43_model()
    model_flux = np.interp(wav, model_wav, model_flux)*(bin_hi-bin_lo)

    expected_rate = np.sum(model_flux / 1.602e-9 / w2e(wav) * specresp)

    print(expected_rate)


if __name__ == '__main__':
    main()
