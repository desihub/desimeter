import numpy as np
from desimeter.posparams.movemask import movemask


def eval_move_flags(posintT, posintP, ptlX, ptlY) :
    '''Simple test to verify if a positioner is moving correctly.

    INPUTS:
        posintT  ... list of theta angles as internally-tracked (a.k.a. POS_T, a.k.a. t_ext)
        posintP  ... list of phi angles as internally-tracked (a.k.a. POS_P, a.k.a. p_ext)

        ptlX     ... list of x as meas with FVC and converted to petal coords (a.k.a. X_PTL)
        ptlY     ... list of y as meas with FVC and converted to petal coords (a.k.a. Y_PTL)

    OUTPUTS:
        flags (composed of movemask bits)
    '''

    flags=0


    assert len(posintT) == len(posintP) == len(ptlX) == len(ptlY)
    assert all(isinstance(arg, list) for arg in (posintT, posintP, ptlX, ptlY))

    t_int = np.array(posintT)
    p_int = np.array(posintP)
    x_ptl = np.array(ptlX)
    y_ptl = np.array(ptlY)


    # loop over unique values of intP to check T moves
    unique_p_int = np.unique(p_int)

    tested = False
    moving = False
    print("DEBUG: unique P values = {}".format(unique_p_int))
    for val in unique_p_int :
        selection = np.where(p_int==val)[0]
        if len(selection)<4 :
            print("DEBUG: for p={} , {} pts: too few points to detect anything".format(val,selection.size))
            continue # too few points to detect anything
        if val>160 :
            print("DEBUG: p={} phi arm is too much retracted to measure variations in theta".format(val))
            continue # phi arm is too much retracted to measure variations in theta
        t_int_range = np.max(t_int[selection])-np.min(t_int[selection])
        if t_int_range<20. :
            print("DEBUG: range of t={}: range of theta angle is to small".format(t_int_range))
            continue # range of theta angle is to small

        tested = True
        rms_xy = np.sqrt(np.std(x_ptl[selection])**2+np.std(y_ptl[selection])**2)
        moving |= (rms_xy>0.1) # this is super arbitrary
    if tested and ( not moving ) :
        print("WARNING: not moving in theta")
        flags |= movemask['THETA_STUCK']
    if not tested :
        print("WARNING: cannot test theta with this dataset")

    return flags
