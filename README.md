# Redshifting

Use this package to shift radio sources in a FITS image from a given redshift to a higher redshift in image plane. \
Be aware that this does also modifies the background RMS.

You can install this package by:
```pip install redshifting```

Example to move the Perseus cluster from redshift 0.018 to 1.318:

```
from redshifting import move_source
myfits = 'Perseus_HBA_full.fits'
move_source(input_fits=myfits, orig_z=0.018, dz=1.3, output_fits='shifted_'+myfits)
```

Video that domstrates the redshifting:
https://youtube.com/shorts/yTT-3RXxbEc
