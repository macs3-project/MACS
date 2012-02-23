available_c = ['0.05','0.01']
row_size = 1001
table_size = row_size*row_size

from array import array as pyarray
import os

class InvalidCError(Exception):
    """Exception about invalid C.

    e.g.:

    raise InvalidCError(0.6)
    """
    def __init__ (self, c):
        self.c = str(c)
        
    def __str__ (self):
        return repr( "Specified C is not available: %s; Should be %s" % (self.strand, "\,".join(available_c)) )

class PreCompiledGFold:
    def __init__ ( self, c ):
        path = os.path.dirname(__file__)
        c = str(c)
        if c in available_c:
            self.gfolds = pyarray('f',[])
            self.gfolds.fromfile(file(os.path.join(path,'g'+str(c)+'.dat'),'rb'),table_size)
        else:
            raise InvalidCError(c)
    
    def get ( self, v1, v2 ):
        i = int(round(v1)*row_size+round(v2))
        return self.gfolds[i]

class FakePreCompiledGFold:
    def __init__ ( self, c ):
        self.gfolds = None
    
    def get ( self, v1, v2 ):
        raise IndexError
    
