#!/usr/bin/env python
# coding: utf-8

import numpy as np
import xarray as xr

from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument("-p", dest="pos_data", required=True)
parser.add_argument("-n", dest="neg_data", required=True)
parser.add_argument('-o', dest='outfile', help='Output file', required=True)


args = parser.parse_args()

# Load the data
#data = xr.open_dataset(args.data)
neg = xr.open_dataset(args.neg_data)
pos = xr.open_dataset(args.pos_data)

data = xr.concat([pos, neg], dim='position')

y = np.concatenate([np.ones (pos.dims['position'], dtype=bool), 
                    np.zeros(neg.dims['position'], dtype=bool)])

data['label'] = xr.DataArray(y, dims = 'position', 
                                coords = dict(position=data.position))

print (f'Writing combined data to {args.outfile}..')
data.to_netcdf(args.outfile)