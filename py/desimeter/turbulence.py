import numpy as np
import astropy.io.ascii
from scipy import optimize
import scipy.linalg
from desimeter.match_positioners import match


def make_covar_gradwavefront(data, param, rq=False):
    sigma, aa, ll = param
    xdist = (data['x'][None, :] - data['x'][:, None])/ll
    ydist = (data['y'][None, :] - data['y'][:, None])/ll
    gauss = aa**2*np.exp(-(xdist**2+ydist**2)/2)
    nn = len(data)
    covar = np.empty((nn*2, nn*2), dtype='f4')
    covar[:nn, :nn] = (1-xdist*xdist)*gauss
    covar[nn:, :nn] = -xdist*ydist*gauss
    covar[:nn, nn:] = -xdist*ydist*gauss
    covar[nn:, nn:] = (1-ydist*ydist)*gauss
    covar += np.diag(np.full(2*nn, sigma**2))
    return covar


def loss_gradwavefront(data, param, rq=False):
    covar = make_covar_gradwavefront(data, param, rq=rq)
    chol, low = scipy.linalg.cho_factor(covar, check_finite=False,
                                        overwrite_a=True)
    dvec = np.concatenate([data['dx'], data['dy']])
    cinvd = scipy.linalg.cho_solve((chol, low), dvec, check_finite=False)
    chi2d = np.dot(dvec, cinvd)
    logdetcovar = 2*np.sum(np.log(np.diag(chol)))
    chi2 = chi2d + logdetcovar
    return 0.5*chi2


def make_covar_independent(data, param, rq=False):
    # might be good to have more large-scale coherence
    # "rational quadratic" seems like the next obvious choice.
    if not rq:
        sigma, aa, ll = param
    else:
        sigma, aa, ll, alpha = param
        alpha = np.clip(alpha, 0.1, 1000)
    dist2 = ((data['x'][None, :] - data['x'][:, None])**2 +
             (data['y'][None, :] - data['y'][:, None])**2)/ll**2
    if not rq:
        covar = aa**2*np.exp(-dist2/2)
    else:
        covar = aa**2*(1+dist2/(2*alpha))**(-alpha)
    covar += np.eye(len(data))*sigma**2
    return covar


def loss_independent(data, param, rq=False):
    covar = make_covar_independent(data, param, rq=rq)
    # chol = np.linalg.cholesky(covar)
    chol, low = scipy.linalg.cho_factor(covar, check_finite=False,
                                        overwrite_a=True)
    cinvx = scipy.linalg.cho_solve((chol, low), data['dx'], check_finite=False)
    cinvy = scipy.linalg.cho_solve((chol, low), data['dy'], check_finite=False)
    chix2 = np.dot(data['dx'], cinvx)
    chiy2 = np.dot(data['dy'], cinvy)

    # factor of 2 for cholesky, and another factor of 2 for x & y
    logdetcovar = 2*2*np.sum(np.log(np.diag(chol)))
    chi2 = chix2+chiy2+logdetcovar
    return 0.5*chi2


def make_data(expect, measure):
    mex, mme = match(expect['LOCATION'], measure['LOCATION'])
    data = np.zeros(len(mex), dtype=[
        ('location', 'i4'), ('x', 'f8'), ('y', 'f8'),
        ('dx', 'f8'), ('dy', 'f8')])
    data['location'] = expect['LOCATION'][mex]
    data['x'] = expect['X_FP'][mex]
    data['y'] = expect['Y_FP'][mex]
    data['dx'] = measure['X_FP'][mme] - data['x']
    data['dy'] = measure['Y_FP'][mme] - data['y']
    return data


def solve_covar(data, lossfun, covarfun, rq=False, nuse=500, **kw):
    if nuse > 0:
        rr2 = data['x']**2+data['y']**2
        s = np.argsort(rr2)
        datause = data[s[:nuse]]
    else:
        datause = data

    def loss(param):
        return lossfun(datause, param, rq=rq)

    # sigma, A, l
    guess = [0.01, 0.01, 10]
    if rq:
        guess = guess + [2]
    res = optimize.minimize(loss, guess, **kw)
    covar = covarfun(data, res.x, rq=rq)
    return covar, res


def solve_files(expectedfn, measuredfn, mode='independent', **kw):
    expect = astropy.io.ascii.read(expectedfn)
    measure = astropy.io.ascii.read(measuredfn)
    data = make_data(expect, measure)
    if mode == 'independent':
        xturb, yturb, res = solve_independent(data, **kw)
    elif mode == 'gradwavefront':
        xturb, yturb, res = solve_gradwavefront(data, **kw)
    else:
        raise ValueError('unknown mode')
    out = measure.copy()
    md, mo = match(data['location'], out['LOCATION'])
    newnames = ['XNEW', 'YNEW', 'DXTURB', 'DYTURB']
    for name in newnames:
        out[name] = np.nan
    out['XNEW'][mo] = out['X_FP'][mo]-xturb[md]
    out['YNEW'][mo] = out['Y_FP'][mo]-yturb[md]
    out['DXTURB'][mo] = xturb[md]
    out['DYTURB'][mo] = yturb[md]
    return out, res


def correct(x, y, x0, y0, dx=None, dy=None):
    data = np.zeros(len(x), dtype=[
        ('x', 'f8'), ('y', 'f8'),
        ('dx', 'f8'), ('dy', 'f8')])
    data['x'] = x
    data['y'] = y
    data['dx'] = x-x0
    data['dy'] = y-y0
    xturb, yturb, _ = solve_independent(data, nuse=500, excludeself=True)
    return x-xturb, y-yturb


def solve_independent(data, excludeself=False, **kw):
    covar, res = solve_covar(data, lossfun=loss_independent,
                             covarfun=make_covar_independent, **kw)

    cninv = np.eye(len(data))*res.x[0]**(-2)
    if not excludeself:
        cpcninv = np.dot(covar, cninv)
        aa = cpcninv+np.eye(len(data))
        xturb = np.linalg.solve(aa, np.dot(cpcninv, data['dx']))
        yturb = np.linalg.solve(aa, np.dot(cpcninv, data['dy']))
    else:
        # Rasmussen & Williams 5.12
        kinv = np.linalg.inv(covar)
        xturb = data['dx']-kinv.dot(data['dx'])/np.diag(kinv)
        yturb = data['dy']-kinv.dot(data['dy'])/np.diag(kinv)
    return xturb, yturb, res


def solve_gradwavefront(data, excludeself=False, **kw):
    covar, res = solve_covar(data, lossfun=loss_gradwavefront,
                             covarfun=make_covar_gradwavefront, **kw)

    dvec = np.concatenate([data['dx'], data['dy']])
    cninv = np.eye(len(dvec))*res.x[0]**(-2)
    if not excludeself:
        cpcninv = np.dot(covar, cninv)
        aa = cpcninv+np.eye(len(dvec))
        turb = np.linalg.solve(aa, np.dot(cpcninv, dvec))
    else:
        # Rasmussen & Williams 5.12
        kinv = np.linalg.inv(covar)
        turb = dvec - kinv.dot(dvec)/np.diag(kinv)
    xturb, yturb = turb[:len(data)], turb[len(data):]
    return xturb, yturb, res


def empirical_covariance(data, bins=10, edges=None):
    dist = np.sqrt((data['x'][None, :] - data['x'][:, None])**2 +
                   (data['y'][None, :] - data['y'][:, None])**2).ravel()
    # meandx = np.mean(data['dx'])
    # meandy = np.mean(data['dy'])
    meandx = 0
    meandy = 0
    dx1dx2 = ((data['dx'][None, :] - meandx) *
              (data['dx'][:, None] - meandx)).ravel()
    dy1dy2 = ((data['dy'][None, :] - meandy) *
              (data['dy'][:, None] - meandy)).ravel()
    if edges is not None:
        bins = len(edges)-1
    else:
        edges = np.linspace(np.min(dist), np.max(dist)+0.1, bins+1)
    ind = np.searchsorted(edges, dist, side='right')-1
    m = ind >= 0
    ncount = np.bincount(ind[m], minlength=bins)
    ncount = ncount + (ncount == 0)
    dx1dx2 = np.bincount(ind[m], weights=dx1dx2[m], minlength=bins)
    dy1dy2 = np.bincount(ind[m], weights=dy1dy2[m], minlength=bins)
    return (edges[:-1]+edges[1:])/2, dx1dx2/ncount, dy1dy2/ncount


def turbulence_gallery(fn, expectfn):
    from matplotlib import pyplot as p
    p.clf()
    expect = astropy.io.ascii.read(expectfn)
    p.subplots_adjust(hspace=0., wspace=0.)
    for i, fn0 in enumerate(fn):
        measure = astropy.io.ascii.read(fn0)
        data = make_data(expect, measure)
        covar, res = solve_covar(data, lossfun=loss_independent,
                                 covarfun=make_covar_independent)
        print(res.x)
        uu, ss, _ = np.linalg.svd(covar)
        xvec = np.dot(uu, np.random.randn(len(ss))*np.sqrt(ss))
        yvec = np.dot(uu, np.random.randn(len(ss))*np.sqrt(ss))
        for k, (dx, dy) in enumerate([
                [data['dx'], data['dy']], [xvec, yvec]]):
            p.subplot(2, 3, 2*i+k+1)
            p.quiver(data['x'], data['y'], dx, dy,
                     units='x', scale=0.001)
            p.gca().xaxis.set_ticklabels('')
            p.gca().yaxis.set_ticklabels('')
            p.gca().xaxis.set_ticks([])
            p.gca().yaxis.set_ticks([])
            p.gca().set_aspect('equal')
