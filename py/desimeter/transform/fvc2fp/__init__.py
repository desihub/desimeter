#- default transform fit
def fit(spots, update_spots=False):
    from .zb import FVCFP_ZhaoBurge
    tx = FVCFP_ZhaoBurge()
    tx.fit(spots, update_spots=update_spots)
    return tx

#- Read JSON file and determine what class to load
def read_jsonfile(filename):
    with open(filename) as fx:
        s = fx.read()

    if s['method'] == 'Zhao-Burge':
        from .zb import FVCFP_ZhaoBurge
        return FVCFP_ZhaoBurge.fromjson(s)
    elif s['method'] == 'xy polynomial':
        from .poly2d import FVCFP_Polynomial
        return FVCFP_Polynomial.fromjson(s)
    else:
        raise ValueError('Unknown method {}'.format(s['method']))
