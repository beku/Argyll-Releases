
This directory has code for taking raw device
information (created from, SPECTRO or CHART), and
extracting an ICC device profile from it.


Profile.exe:
------------

profile.exe	takes the target patch values, and creates
and ICC profile.

Options: 

 -v              Verbose mode

 -q [lmh]        Quality - Low, Medium (def), High

 -b				 Low quality B2A table

 -y              Verify A2B profile

 -k [hzp]        h = 0.5 k (default)
                 z = zero k

 -k p black white npow sat spow
                 p = shape parameters

 -l limit        set ink limit, 1 - 400%

 -d [lgs]        Display profile type
                 l = lut, g = gamma+matrix, s = shaper+matrix

 -i illum        Choose illuminant for print/transparency spectral data:
                 D50 (def.), D65, F8, F10

 -o observ       Choose CIE Observer for spectral data:
                 1931_2 (def.), 1964_10, 1978_2

 outfile         Base name for input.ti3/output.icm file


For CMYK or RGB output devices, only Lut based profiles can be created.

For RGB Display devices, matrix based profiles can be created.
