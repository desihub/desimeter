import numpy as np
import astropy.io.ascii
from scipy import optimize
import scipy.linalg
from desimeter.match_positioners import match
import scipy.spatial


# accommodate both use in online system and in datasystems environment
# use online system logging when available
try:
    import DOSLib.logger as log
    log.warning = log.warn
except Exception:
    from desimeter.log import get_logger
    log = get_logger()


# =================================================================
# from Sergey : correction using local polynomial fit
# =================================================================

def getpoly(x, y, ndeg=2):
    """
    Get the 2D polynomial design matrix
    """
    N = len(x)
    polys = np.zeros((((ndeg + 1) * (ndeg + 2)) // 2 - 1, N * 2))  # []
    cnt = 0
    for deg in range(1, ndeg + 1):
        for j in range(deg + 1):
            i1, i2 = j, deg - j
            # notice the i1==0 i2==0 are there to prevent 0**(-1)
            # these are dF/dx dF/dy of the polynomial
            polys[cnt, :N] = i1 * x**(i1 - 1 + (i1 == 0)) * y**i2
            polys[cnt, N:] = i2 * x**i1 * y**(i2 - 1 + (i2 == 0))
            cnt += 1
    return polys


def predictor(curx, cury, curdx, curdy, tofit, ndeg=2, polys=None):
    """
    Fit the offsets by a polynomial f-n
    return the model predictions and
    the cross-validated norm
    """
    ncv = 3
    if polys is None:
        polys = getpoly(curx, cury, ndeg=ndeg)
    else:
        polys = polys[:((ndeg + 1) * (ndeg + 2)) // 2 - 1]
    cvid = np.arange(2 * len(curx)) % ncv
    norm = 0
    dxy = np.concatenate((curdx, curdy))
    tofit2 = np.concatenate((tofit, tofit))
    for i in range(ncv):
        cursub = tofit2 & (cvid != i)
        cursub1 = tofit2 & (cvid == i)
        xcoeff = scipy.linalg.basic.lstsq(polys[:, cursub].T,
                                          dxy[cursub],
                                          check_finite=False)[0]
        xpred = np.dot(polys.T, xcoeff)[cursub1]
        norm += np.sum((xpred - dxy[cursub1])**2)
    xcoeff = scipy.linalg.basic.lstsq(polys[:, tofit2].T, dxy[tofit2])[0]
    xpred, ypred = np.dot(polys.T, xcoeff)[~tofit2]
    return xpred, ypred, norm, polys


def correct_with_pol(x, y, x0, y0, win=50):
    """
    Parameters
    ----------
    x: ndarray
        Measured x
    y: ndarray
        Measured y
    x0: ndarray
        Expected x
    y0: ndarray
        Expected y

    Returns
    -------
    xy: tuple of ndarray
        Tuple of corrected arrays

    """
    maxndeg = 4
    # offsets wrt reference
    dx = x - x0
    dy = y - y0

    X0 = np.array([x0, y0]).T
    X = np.array([x, y]).T

    T0 = scipy.spatial.cKDTree(X0)
    N = len(x)

    # predicted offset
    dxpred = np.zeros(N)
    dypred = np.zeros(N)
    bestdegs = np.zeros(N)

    for i in range(N):
        # go over each point
        # query the neighborhood
        xids = T0.query_ball_point(X[i], win)
        xids = np.array(xids)

        curxcen, curycen = x[i], y[i]
        curx = x[xids] - curxcen
        cury = y[xids] - curycen
        tofit = xids != i
        curdx = dx[xids]
        curdy = dy[xids]
        bestnorm = 1e9
        # try polynomials of different degrees
        # selecting by cross-validate norm
        polys = None
        for ndeg in range(maxndeg, 0, -1):
            xpred, ypred, norm, polys = predictor(curx,
                                                  cury,
                                                  curdx,
                                                  curdy,
                                                  tofit,
                                                  ndeg=ndeg,
                                                  polys=polys)
            if norm < bestnorm:
                bestnorm = norm
                lastx, lasty = xpred, ypred
                bestdeg = ndeg
        dxpred[i] = lastx
        dypred[i] = lasty
        bestdegs[i] = bestdeg
    return x - dxpred, y - dypred
# ===============================================================


# ===============================================================
# From Eddie : correction using gaussian processes
# ===============================================================


def make_covar_gradwavefront_nonoise(x1, y1, x2, y2, param, rq=False, **kw):
    if rq:
        raise ValueError('gradwavefront does not support rq mode.')
    _, aa, ll = param
    xdist = (x1[None, :] - x2[:, None])/ll
    ydist = (y1[None, :] - y2[:, None])/ll
    gauss = aa**2*np.exp(-(xdist**2+ydist**2)/2)
    n1 = len(x1)
    n2 = len(x2)
    covar = np.empty((n2*2, n1*2), dtype='f4')
    covar[:n2, :n1] = (1-xdist*xdist)*gauss
    covar[n2:, :n1] = -xdist*ydist*gauss
    covar[:n2, n1:] = -xdist*ydist*gauss
    covar[n2:, n1:] = (1-ydist*ydist)*gauss
    return covar


def make_covar_gradwavefront(data, param, rq=False, **kw):
    """Make derivative-of-Gaussian covariance matrix for data given param.

    Args:
        data : ndarray containing x, y, dx, dy for each positioner
        param : array_like(float) giving [sigma, amplitude, length scale]
            for Gaussian covariance

    Returns:
        covar : covariance matrix (n*len(data)*2, n*len(data)*2) for [dx, dy]
    """
    covar = make_covar_gradwavefront_nonoise(
        data['x'], data['y'], data['x'], data['y'], param, rq=rq)
    sigma = param[0]
    covar += np.diag(np.full(covar.shape[0], sigma**2))
    return covar


def loss_gradwavefront(data, param, rq=False):
    """Loss function for data using derivative-of-Gaussian covariance matrix.

    Args:
        data : ndarray containing x, y, dx, dy for each positioner
        param : array_like(float) giving [sigma, amplitude, length scale]
            for Gaussian covariance

    Returns:
        loss : -0.5 * log(marginal likelihood) of data coming from the given
           covariance matrix
    """
    covar = make_covar_gradwavefront(data, param, rq=rq)
    chol, low = scipy.linalg.cho_factor(covar, check_finite=False,
                                        overwrite_a=True)
    dvec = np.concatenate([data['dx'], data['dy']])
    cinvd = scipy.linalg.cho_solve((chol, low), dvec, check_finite=False)
    chi2d = np.dot(dvec, cinvd)
    logdetcovar = 2*np.sum(np.log(np.diag(chol)))
    chi2 = chi2d + logdetcovar
    return 0.5*chi2


def make_covar_independent_nonoise(x1, y1, x2, y2, param, rq=False, **kw):
    if not rq:
        _, aa, ll = param
    else:
        _, aa, ll, alpha = param
        alpha = np.clip(alpha, 0.1, 1000)
    dist2 = ((x1[None, :] - x2[:, None])**2 +
             (y1[None, :] - y2[:, None])**2)/ll**2
    if not rq:
        covar = aa**2*np.exp(-dist2/2)
    else:
        covar = aa**2*(1+dist2/(2*alpha))**(-alpha)
    return covar


def make_covar_independent(data, param, rq=False, **kw):
    """Gaussian covariance matrix for data given param.

    Args:
        data : ndarray containing x, y, dx, dy for each positioner
        param : array_like(float) giving [sigma, amplitude, length scale]
            for Gaussian covariance
        rq : bool, use rational quadratic covariance instead of Gaussian

    Returns:
        covar : covariance matrix for dx (equivalently, dy)
    """
    covar = make_covar_independent_nonoise(
        data['x'], data['y'], data['x'], data['y'], param, rq=rq)
    sigma = param[0]
    covar += np.eye(len(data))*sigma**2
    return covar


def loss_independent(data, param, rq=False):
    """Loss function for data given param describing covariance matrix.

    Args:
        data : ndarray containing x, y, dx, dy for each positioner
        param : array_like(float) giving [sigma, amplitude, length scale]
            for Gaussian covariance

    Returns:
        loss : -0.5 * log(marginal likelihood) of data coming from the given
           covariance matrix
    """
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
    """Convert expected and measured desimeter data structures to data array.

    Args:
        expect : desimeter data structure giving expected locations of
            positioners.
        measure : desimeter data structure giving measured locations of
            positioners.

    Returns:
        ndarray containing (location, x, y, dx, dy) for positioners in
            expect and measure.
    """
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
    """Find covariance matrix best explaining data.

    Args:
        lossfun : function taking data & param, returning -0.5 log(likelihood)
        covarfun : function taking data & param, returning covariance matrix
        rq : bool, use rational quadratic covariance matrix rather than
            Gaussian
        nuse : int, use only central nuse fibers when optimizing covariance
            matrix
        **kwargs : additional arguments passed to scipy.optimize.minimize

    Returns:
        covar : ndarray, best fit covariance matrix
        res : output of scipy.optimize.minimize
    """
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
    """Find turbulent contribution to measured positions of fibers.

    Args:
        expectedfn : file name containing expected locations of positioners
        measuredfn : file name containing measured locations of positioners
        mode : 'independent' or 'gradwavefront'
            models turbulent offset correlations as either a Gaussian process
            with a Gaussian correlation function or a derivative-of-Gaussian
            correlation function (appropriate if the turbulent contributions
            can be modeled as the derivative of a Gaussian-correlated wavefront
            offset)
        **kw : additional keywords passed to solve_{independent, gradwavefront}

    Returns:
        loss : -0.5 * log(marginal likelihood) of data coming from the given
           covariance matrix
    """
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


def correct_using_stationary(xs, ys, x0s, y0s, xc, yc,
                             scale_covar=1):
    """Remove correlated turbulent signals from measured locations of fibers,
    using subset of stationary fibers.

    See `correct` for more documentation.  This function extends `correct`
    to use a subset of the fibers to correct all the fibers.

    Args:
        xs : array_like(n), measured x position of stationary fibers
        ys : array_like(n), measured y position of stationary fibers
        x0s : array_like(n), expected x position of stationary fibers
        y0s : array_like(n), expected y position of stationary fibers
        xc : array_like(n), measured x position of fibers to be corrected
        yc : array_like(n), measured y position of fibers to be corrected
        scale_covar : bool, amount to scale covariance by
            (pix / micron for FVC instead of FP, e.g.)

    Returns:
        xn, yn: turbulence-corrected positions of fibers
    """
    data = np.zeros(len(xs), dtype=[
        ('x', 'f8'), ('y', 'f8'),
        ('dx', 'f8'), ('dy', 'f8')])
    data['x'] = xs
    data['y'] = ys
    data['dx'] = xs-x0s
    data['dy'] = ys-y0s
    xturb, yturb, _ = solve_independent(
        data, nuse=500, excludeself=False,
        predict_at=(xc, yc), method='powell',
        fix_covar=True, scale_covar=scale_covar)
    return xc-xturb, yc-yturb


def correct(x, y, x0, y0, dx=None, dy=None):
    """Remove correlated turbulent signals from measured locations of fibers.

    Residuals x-x0, y-y0 are modeled as a Gaussian process.  The
    covariance of this Gaussian process is fit as the sum of the
    following components:
    - a diagonal measurement error and positioner movement error term
    - a spatially correlated turbulent term

    The spatially correlated turbulent term is modeled with a Gaussian
    covariance matrix whose parameters are fit using the inner 500 positioners.
    The covariance is assumed to be isotropic in x & y, but the best fit
    turbulence is applied independently to x & y.

    dx & dy are placeholders for future code that includes different
    uncertainties for different fibers.

    Args:
        x : array_like(n), measured x positions of fibers
        y : array_like(n), measured y positions of fibers
        x0 : array_like(n), expected x positions of fibers
        y0 : array_like(n), expected y positions of fibers
        dx: array_like(n), uncertainty in x.  Not currently used.
        dy: array_like(n), uncertainty in y.  Not currently used.

    Returns:
        x : array_like(n), turbulence-corrected x positions of fibers
        y : array_like(n), turbulence-corrected y positions of fibers
    """
    if (dx is not None) or (dy is not None):
        log.warning('Uncertainties in x & y are currently ignored.')
    data = np.zeros(len(x), dtype=[
        ('x', 'f8'), ('y', 'f8'),
        ('dx', 'f8'), ('dy', 'f8')])
    data['x'] = x
    data['y'] = y
    data['dx'] = x-x0
    data['dy'] = y-y0
    xturb, yturb, _ = solve_independent(data, nuse=500, excludeself=True)
    return x-xturb, y-yturb


def solve_independent(data, excludeself=False, predict_at=None,
                      fix_covar=False, scale_covar=1, **kw):
    """Find turbulent contributions to measured fiber positions.

    Assumes that the covariance matrix should be the same for x & y,
    but that otherwise the x & y directions are independent.

    Args:
        data : ndarray containing measured positions and residuals
            from expected locations
        excludeself : bool
            do not use this fiber when computing the turbulence
            affecting this fiber.
        predict_at : tuple(array_like, array_like)
            x & y coordinates of locations at which to predict turbulence
        **kw : additional keywords passed to solve_covar

    Returns:
        xturb : turbulent contributions in x direction
        yturb : turbulent contributions in y direction
        res : output from scipy.optimize.minimize describing best fit
            covariance matrix
    """
    if not fix_covar:
        covar, res = solve_covar(data, lossfun=loss_independent,
                                 covarfun=make_covar_independent, **kw)
    else:
        from types import SimpleNamespace
        res = SimpleNamespace()
        res.x = [5e-3*scale_covar, 5e-3*scale_covar, 50*scale_covar]
        if kw.get('rq', False):
            res.x = res.x + [2]
        covar = make_covar_independent(data, res.x, **kw)

    if predict_at is not None and excludeself:
        raise ValueError('predict_at does not make sense in combination with '
                         'excludeself')

    if not excludeself:
        if predict_at:
            # K(X*, X)(K(X, X) + C_n)^-1 y
            # Rasmussen & Williams algorithm 2.1
            chol, low = scipy.linalg.cho_factor(covar, check_finite=False,
                                                overwrite_a=True)
            covarpred = make_covar_independent_nonoise(
                data['x'], data['y'], predict_at[0], predict_at[1],
                res.x, **kw)
            alphax = scipy.linalg.cho_solve((chol, low), data['dx'])
            alphay = scipy.linalg.cho_solve((chol, low), data['dy'])
            xturb = np.dot(covarpred, alphax)
            yturb = np.dot(covarpred, alphay)
        else:
            # covar includes measurement uncertainties; separate that part out.
            cninv = np.eye(len(data))*res.x[0]**(-2)
            covar -= np.eye(len(data))*res.x[0]**2
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


def solve_gradwavefront(data, excludeself=False, predict_at=None,
                        fix_covar=False, **kw):
    """Find turbulent contributions to measured fiber positions.

    Assumes that the turbulent contributions can be modeled as the
    gradient of a wavefront error.  i.e., they are curl free.

    Args:
        data : ndarray containing measured positions and residuals
            from expected locations
        excludeself : bool
            do not use this fiber when computing the turbulence
            affecting this fiber.
        **kw : additional keywords passed to solve_covar

    Returns:
        xturb : turbulent contributions in x direction
        yturb : turbulent contributions in y direction
        res : output from scipy.optimize.minimize describing best fit
            covariance matrix
    """

    if predict_at is not None and excludeself:
        raise ValueError('predict_at does not make sense in combination with '
                         'excludeself')

    if not fix_covar:
        covar, res = solve_covar(data, lossfun=loss_gradwavefront,
                                 covarfun=make_covar_gradwavefront, **kw)
    else:
        from types import SimpleNamespace
        res = SimpleNamespace()
        res.x = [5e-3, 5e-3, 100]
        if kw.get('rq', False):
            res.x = res.x + [2]
        covar = make_covar_gradwavefront(data, res.x, **kw)

    dvec = np.concatenate([data['dx'], data['dy']])
    if not excludeself:
        if predict_at:
            # K(X*, X)(K(X, X) + C_n)^-1 y
            # Rasmussen & Williams algorithm 2.1
            chol, low = scipy.linalg.cho_factor(covar, check_finite=False,
                                                overwrite_a=True)
            covarpred = make_covar_gradwavefront_nonoise(
                data['x'], data['y'], predict_at[0], predict_at[1],
                res.x, **kw)
            alpha = scipy.linalg.cho_solve((chol, low), dvec)
            turb = np.dot(covarpred, alpha)
            xturb, yturb = turb[:len(predict_at[0])], turb[len(predict_at[0]):]
        else:
            # remove measurement noise contribution to covar
            cninv = np.eye(len(dvec))*res.x[0]**(-2)
            covar -= np.eye(len(dvec))*res.x[0]**2
            cpcninv = np.dot(covar, cninv)
            aa = cpcninv+np.eye(len(dvec))
            turb = np.linalg.solve(aa, np.dot(cpcninv, dvec))
            xturb, yturb = turb[:len(data)], turb[len(data):]
    else:
        # Rasmussen & Williams 5.12
        kinv = np.linalg.inv(covar)
        turb = dvec - kinv.dot(dvec)/np.diag(kinv)
        xturb, yturb = turb[:len(data)], turb[len(data):]
    return xturb, yturb, res


def empirical_covariance(data, bins=10, edges=None):
    """Compute empirical covariance of the data."""
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
    """Make some plots of the turbulence in measured data from desimeter."""
    from matplotlib import pyplot as p
    p.clf()
    expect = astropy.io.ascii.read(expectfn)
    p.subplots_adjust(hspace=0., wspace=0.)
    for i, fn0 in enumerate(fn):
        measure = astropy.io.ascii.read(fn0)
        data = make_data(expect, measure)
        covar, _ = solve_covar(data, lossfun=loss_independent,
                                 covarfun=make_covar_independent)
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
