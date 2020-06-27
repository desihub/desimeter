import numpy as np
import desimeter.log

log = desimeter.log.get_logger()

def match2d(x1, y1, x2, y2, rad):
    """Find all matches between x1, y1 and x2, y2 within radius rad."""
    from scipy.spatial import cKDTree
    xx1 = np.stack([x1, y1], axis=1)
    xx2 = np.stack([x2, y2], axis=1)
    tree1 = cKDTree(xx1)
    tree2 = cKDTree(xx2)
    res = tree1.query_ball_tree(tree2, rad)
    lens = [len(r) for r in res]
    m1 = np.repeat(np.arange(len(x1), dtype='i4'), lens)
    if sum([len(r) for r in res]) == 0:
        m2 = m1.copy()
    else:
        m2 = np.concatenate([r for r in res if len(r) > 0])
    d12 = np.sqrt(np.sum((xx1[m1, :]-xx2[m2, :])**2, axis=1))
    return m1, m2, d12


def match(a, b):
    """Find indices ma, mb such that a[ma] == b[mb].  a must have no duplicates.
    """
    sa = np.argsort(a)
    if len(np.unique(a)) != len(a):
        raise ValueError('All keys in a must be unique.')
    ind = np.searchsorted(a[sa], b)
    m = (ind >= 0) & (ind < len(a))
    matches = a[sa[ind[m]]] == b[m]
    m[m] &= matches
    return sa[ind[m]], np.flatnonzero(m)


def match_positioners(fvc, metr, calib, countthresh=80000,
                      arm_fudge_offset=0.1, return_alternatives=False):
    """Take FVC FP measurements, metrology, and calibration measurements,
    figure out best solution.

    To use Stephen's terminology, this plays the ``sudoku'' game ~optimally,
    finding the best possible assignment of positioners to centroids.  It then
    finds a number of equally good alternatives, if they exist.  It returns
    the best assignment for each positioner, a flag indicating if that assignment
    was unique, and optionally the equally good alternatives that were found.

    There is no guarantee that all equally good alternatives were found, only
    that if a positioner had an ambiguous assignment, at least one alternative
    was found where that positioner had an alternative assignment.

    Args:
    ----
        fvc : desimeter fvc structure with info on centroids
        metr : desimeter metrology structure with info on positioners
        calib : desimeter calibration strucuter with info on arm lengths
        countthresh : only fvc centroids with flux > countthresh are fit
        arm_fudge_offset : make all arms longer by arm_fudge_offset.
            This is intended to compensate for uncertainties in the measurements.

    Returns
    -------
    ind, score, alternatives
    ind : array of length fvc, giving the indices into metr of the positioners
        assigned to each centroid in fvc.
    score : the score of this assignment.
    alternatives : if alternatives is set, an array [n_alternatives, n_fvc]
      of alternatives that are equally good assignments.


    Notes
    -----
    This code finds optimal solutions in the limit that there is no uncertainty
    in the arm length and that arm_fudge_offset is 0, and that the only
    geometric constraint on the location of positioners is that they can't reach
    things outside of their arm lengths.  Future improvements to optimality
    could think harder about the actual shapes of the positioners (presumably
    some assignments this code makes are impossible because the positioners
    would need to collide to achieve them), and better incoporate the
    uncertainties in R1+R2.  Another good feature would be to add expected
    positioner/fiber locations, and increase the score of an assignment if
    - there is a centroid in the expected location
    - the correct positioner is assigned to that location
    this could be easily accommodated in the score_matrix below.

    Optimality would also be improved if knowledge like which fibers are
    known to be broken and better exclusions of spurious "centroids" were added.
    Particularly the case of positioners with no arm length measurements
    (presumably broken) could be improved.
    """

    fvcpos = np.flatnonzero((fvc['PINHOLE_ID'] == 0) &
                            (fvc['COUNTS'] > countthresh))
    # need to apply a count threshhold to get rid of some bogus 'centroids'

    metrpos = np.flatnonzero(metr['PINHOLE_ID'] == 0)
    fvclen = len(fvc)
    fvc = fvc[fvcpos]
    metr = metr[metrpos]
    mm, mc = match(metr['DEVICE_ID'], calib['POS_ID'])
    metrarmlen = np.zeros(len(metr), dtype='f4')
    metrarmlen[mm] = (calib['LENGTH_R1_STATIC'][mc] +
                      calib['LENGTH_R2_STATIC'][mc] +
                      arm_fudge_offset)
    # 0.1 is a made up number intended to be higher that most of the
    # uncertainties in R1+R2.  Could be replaced by actual propagation of the
    # uncertainties in R1+R2 multiplied by something like 5.
    # Though I started with 0.06 here, and found cases where that was clearly
    # not enough, even where R1 and R2 had uncertainties of ~0.01.

    # some positioners don't have arm lengths.  I inspected petal the test
    # petal and it looked to me like these were all close to homed.  So I'm
    # artificially going to set their arm lengths to 4.3, short enough that
    # the positioners won't overlap their neighbors, but they can find
    # the homed centroid near them.
    metrarmlen[metrarmlen == 0] = 4.3

    # we have a bunch of positioner center locations in metr and the reach of
    # each positioner.  Now let's find which spots are in reach of which
    # positioners.

    if len(metr) != len(fvc):
        log.info('number of positioners does not equal number of centroids; '
                 'perfect assignment not possible.')
    mm, mf, dmf = match2d(metr['X_FP'].data, metr['Y_FP'].data,
                          fvc['X_FP'].data, fvc['Y_FP'].data,
                          np.max(metrarmlen)+0.1)
    okdistance = np.array([d0 < a0
                           for (d0, a0) in zip(dmf, metrarmlen[mm])])
    mm = mm[okdistance]
    mf = mf[okdistance]
    dmf = dmf[okdistance]

    # okay.  Now we have a all possible matches between positioners and centroids
    # modulo only the uncertainty in the centroid, which is small.
    # what possible assignments of fibers to centroids exist?

    score_matrix = np.zeros((len(metr), len(fvc)), dtype='f4')
    for i, j in zip(mm, mf):
        score_matrix[i, j] = 1

    assignment = solve_assignment(score_matrix)
    nassign = len(assignment[0])

    # we now have the best assignment and score.  Are there other equally
    # good ones?
    # I think this just requires checking the score matrix for cycles, which
    # is very cheap, but I don't want to think about writing that code right now.
    # The other thing we can do is resolve the after removing each edge in the
    # graph.  This is a much more expensive solution, and takes a few minutes
    # for one petal.
    # Note scipy 1.4 is ~1000x faster, and maybe there's time to save if I
    # filter out the unambiguous positioners before this step.

    if return_alternatives:
        alternatives = list()
        new_score_matrix = score_matrix.copy()
        for i in range(nassign):
            if (i % 50) == 0:
                print(f'checking alternative {i} of {nassign}')
            old = new_score_matrix[assignment[0][i], assignment[1][i]]
            new_score_matrix[assignment[0][i], assignment[1][i]] = 0
            alternative = solve_assignment(new_score_matrix)
            alternative[0][:] = metrpos[alternative[0]]
            alternative[1][:] = fvcpos[alternative[1]]
            if alternative[2] == assignment[2]:
                alternatives.append(alternative)
            new_score_matrix[assignment[0][i], assignment[1][i]] = old
        assignment[0][:] = metrpos[assignment[0]]
        assignment[1][:] = fvcpos[assignment[1]]

    ret = np.full(fvclen, -1)
    ret[assignment[1]] = assignment[0]
    res = (ret, assignment[2])

    log.info(f'made {assignment[2]} assignments, between {len(fvc)} '
             f'valid centroids and {len(metr)} valid positioners')

    if return_alternatives:
        ret = np.full((len(alternatives), fvclen), -1)
        for i, alt in enumerate(alternatives):
            ret[i, alt[1]] = alt[0]
        ret = res + (ret,)

    return ret

def solve_assignment(score_matrix):
    from scipy.optimize import linear_sum_assignment
    assignment = linear_sum_assignment(-score_matrix)
    score = np.sum([score_matrix[i, j] for (i, j) in zip(*assignment)])
    return assignment + (score,)


def plot_match(fvc, metr, assignment, alternatives=None):
    possible_assignments = dict()
    if alternatives is not None:
        allassignments = np.vstack([assignment, alternatives])
    else:
        allassignments = assignment[None, ...]
    for a in allassignments:
        for fvcind, metrind in enumerate(a):
            if metrind == -1:
                continue
            possible_assignments[fvcind] = (
                possible_assignments.get(fvcind, set()) | set([metrind]))
    from matplotlib import pyplot as p
    p.plot(fvc['X_FP'], fvc['Y_FP'], '+', label='centroids')
    p.plot(metr['X_FP'], metr['Y_FP'], 'x', label='positioner centers')
    for fvcind in possible_assignments:
        for metrind in possible_assignments[fvcind]:
            p.plot([fvc['X_FP'][fvcind], metr['X_FP'][metrind]],
                   [fvc['Y_FP'][fvcind], metr['Y_FP'][metrind]], 'r-')
    p.gca().set_aspect('equal')
    p.xlabel('X_FP (mm)')
    p.ylabel('X_FP (mm)')
    p.legend()
