from astropy.io import fits

file = "HVC_219_226_-1_9.fits"
hdul = fits.open(file)
hdr = hdul[0].header
data = hdul[0].data

hdr['CRPIX2'] -= hdr['CRVAL2'] / hdr['CDELT2']
hdr['CRVAL2'] = 0
hdr['BUNIT'] = 'K'
hdr['RESTFRQ'] = 1420405751.77
hdr['CTYPE3'] = 'VRAD'

fits.writeto('HVC_219_226_-1_9.fits', data, hdr, overwrite=True)
