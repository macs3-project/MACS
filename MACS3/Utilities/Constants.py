MACS_VERSION = "3.0.0b2"
MAX_PAIRNUM = 1000
MAX_LAMBDA  = 100000
FESTEP      = 20
BUFFER_SIZE = 100000                   # np array will increase at step of 1 million items
READ_BUFFER_SIZE = 10000000            # 10M bytes for read buffer size
N_MP = 2                               # Number of processers

#Effective genome size, collected from
#https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html

EFFECTIVEGS = {"hs":2913022398,                  #GRCh38
               "mm":2652783500,                  #GRCm38
               "ce":100286401,                   #WBcel235
               "dm":142573017}                   #dm6
