# Redshifting

Use this package to shift radio sources in image plane from a given redshift to a higher redshift using FITS images. \
Be aware that this also scales the background RMS.

You can install this package by:
```pip install redshifting```

Example for shifting from redshift 0.018 to 1.318:

```
from redshifting import move_source
myfits = 'Perseus_HBA_full.fits'
move_source(input_fits=myfits, z=0.018, dz=1.3, output_fits='shifted_'+myfits)
```

Video that domstrates the redshifting:
https://youtube.com/shorts/yTT-3RXxbEc

### Citation

If you use this code, please cite:

> de Jong, J.M.G.H.J. et al. (2024), *Cosmic evolution of FRI and FRII sources out to z = 2.5*,  
> Astronomy & Astrophysics, 683, A23.  
> DOI: https://doi.org/10.1051/0004-6361/202347131  
> arXiv: https://arxiv.org/abs/2311.13427