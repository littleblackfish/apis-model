import json
from itertools import product
from random import sample

n = 25

conv_nfilters = [64,128,256]
conv_ksize = [8,16,32]
conv_psize = [8,16,32]
dense_units = [32,64,128]
dropout = [0,0.25,0.5]

pool = list(product(conv_nfilters, conv_ksize, conv_psize, dense_units, dropout))


models = []
for i, params  in enumerate(sample(pool, 25)):
    model = dict(name = f'onehot_convo_{i}',
                 conv_nfilters=params[0],
                 conv_ksize = params[1],
                 conv_psize = params[2],
                 dense_units = params[3],
                 dropout = params[4]
                 )

    models.append(model)

with open('models.json', 'w') as f :
    json.dump(models, f)

