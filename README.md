# SaclayMocks

Simulated Absorption for Cosmology with Lyman-Alpha from the Yvette Mocks


## Installation:
to download, you can simply use git clone:
```
git clone https://github.com/igmhub/SaclayMocks.git
```

then, you need to add to your bashrc:
```
export SACLAYMOCKS_BASE=<path to the git repository>
export PATH=$SACLAYMOCKS_BASE/bin/:$PATH
export PYTHONPATH=$SACLAYMOCKS_BASE/py/:$PYTHONPATH
```

## Running:
The code runs as follow:
First it reads and interpolate the powerspectrum from camb (stored in `etc/`) to create a 3D power-spectrum.
These steps are done with `interpolate_pk.py` and `merge_pk.py`

Then, the Gaussian random field boxes are drawn, with `make_boxes.py`

The quasars are drawn from the GRF boxes, using `draw_qso.py`

Then, the spectra are computed for each quasar, using `make_spectra.py` and `merge_spectra.py`

All theses steps are done for each of the 7 chunks that make up the desi footprint

Once all the chunks are done, they are combined to produce the final output. This is done by running:
`merge_qso.py` for producing the `master.fits` and `master_randoms.fits` files
`dla_saclay.py` and `merge_dla.py` to produce the `master_DLA.fits` file
`dla_randoms.py` to produce the `master_DLA_randoms.fits` file
`make_transmissions.py` to produce all the transmissions files


You don't need to run all these codes by hand. You just have to run `submit_mocks.py`
This code will produce several bash script that will run all the needed python code, according to the options.

To produce one realisation, just run:
```
submit_mocks.py --mock-dir <out_path>
```

and it will produce a bash script `submit.sh`, that will send all the bash script which
themselves send the python scripts to the cori nodes, on NERSC.


If you want to produce a small footprint, on your laptop, you can run:
```
submit_mocks.py --mock-dir <out_path> --cori-nodes False --box-size 256
```
When not running on cori nodes, it is advised to not use the nominal chunk size (2560),
because it uses a lot of memory and threads.
Also, if not distributed on several nodes, it can take hours.

You can also limit the number of chunks by specifying the IDs of chunks to runs:
```
--chunk-id 1 2
```

All the informations relative to chunks (box size available, chunk ids, ...) are in the `chunk_parameters()` function.
