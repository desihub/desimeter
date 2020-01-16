"""
Base class for transformations between FVC and FP coordinates
"""

class FVC2FP_Base(object):
    """
    Base class for transforms between FVC and FP coordinates

    Subclasses implement this interface with specific parameterizations
    """
    def __init__(self):
        pass
    
    def tojson(self):
        raise NotImplementedError
    
    @classmethod
    def fromjson(cls, jsonstring):
        raise NotImplementedError
    
    @classmethod
    def read_jsonfile(cls, filename):
        with open(filename) as fx:
            s = fx.read()
        return cls.fromjson(s)

    def write_jsonfile(self, filename):
        with open(filename, 'w') as fx:
            fx.write(self.tojson())

    def fit(self, spots, metrology=None, update_spots=False):
        raise NotImplementedError

    def fvc2fp(self, xpix, ypix, xerr=None, yerr=None):
        """
        Converts fiber view camera pixel x,y -> focal plane x,y
        """
        raise NotImplmentedError

    def fp2fvc(self, xfp, yfp):
        """
        Converts focal plane x,y -> fiber view camera pixel x,y
        """
        raise NotImplmentedError
