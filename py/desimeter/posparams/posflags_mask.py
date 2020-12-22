'''These are flags set during online operation of the instrument. Referred to as
"pos_flags" within the PetalApp / petal context.

The values here were retrieved by J.Silber on 2020-12-22.
https://desi.lbl.gov/svn/code/online/DOSlib/trunk/python/DOSlib/flags.py
SVN r134618
'''

from desimeter.bitmask import BitMask

POSITIONER_FLAGS_BITS = {'MATCHED' : 0,
                         'PINHOLE' : 1,
                         'POSITIONER' : 2,
                         'FIDUCIAL' : 3,
                         'FVCERROR' : 4,
                         'BADPOSFID' : 5,
                         'MOVED' : 6,
                         'GIF' : 7,
                         'ETC' : 8,
                         'FITTEDPINHOLE' : 9,
                         'MATCHEDCENTER' : 10,
                         'AMBIGUOUS' : 11,
                         'FVCPROC' : 12,
                         'ASSIGNED' : 14,
                         'CONVERGED' : 15,
                         'NOTCTLENABLED' : 16,
                         'BROKENFIBER' : 17,
                         'COMERROR' : 18,
                         'OVERLAP' : 19,
                         'FROZEN' : 20,
                         'UNREACHABLE' : 21,
                         'BOUNDARYVIOLATION' : 22,
                         'MULTIPLEREQUESTS' : 23,
                         'NONFUNCTIONAL' : 24,
                         'REJECTED' : 25,
                         'EXPERTLIMIT' : 26,
                         'BADNEIGHBOR' : 27,
                         'MISSINGSPOT' : 28,
                         'BADPERFORMANCE' : 29
                    }

POSITIONER_FLAGS_VERBOSE = {'MATCHED' : 'Matched',
                            'PINHOLE' : 'Fiducial pinhole',
                            'POSITIONER' : 'Positioner',
                            'FIDUCIAL' : 'Fiducial',
                            'FVCERROR' : 'FVC error',
                            'BADPOSFID' : 'FVC bad posfid',
                            'MOVED' : 'Moved',
                            'GIF' : 'GIF',
                            'ETC' : 'ETC',
                            'FITTEDPINHOLE' : 'Fitted pinhole',
                            'MATCHEDCENTER' : 'Matched to center',
                            'AMBIGUOUS' : 'Ambiguous match',
                            'FVCPROC' : 'FVCPROC',
                            'ASSIGNED' : 'Target assiged',
                            'CONVERGED' : 'Converged on target',
                            'NOTCTLENABLED' : 'Control disabled',
                            'BROKENFIBER' : 'Fiber not intact',
                            'COMERROR' : 'CAN communication error',
                            'OVERLAP' : 'Overlapping targets',
                            'FROZEN' : 'Frozen by anticollision',
                            'UNREACHABLE' : 'Unreachable target',
                            'BOUNDARYVIOLATION' : 'Exceeds boundaries',
                            'MULTIPLEREQUESTS' : 'Multiple requests',
                            'NONFUNCTIONAL' : 'Device nonfunctional',
                            'REJECTED' : 'Movetable rejected',
                            'EXPERTLIMIT' : 'Exceeds expert limit',
                            'BADNEIGHBOR' : 'Bad neighbor',
                            'MISSINGSPOT' : 'Unmatched',
                            'BADPERFORMANCE' : 'Exceeded max error'
                           }

_bitdefs = {'posflags_mask': [[key, POSITIONER_FLAGS_BITS[key], POSITIONER_FLAGS_VERBOSE[key]] for key in POSITIONER_FLAGS_BITS]}

posflags_mask = BitMask('posflags_mask', _bitdefs)
