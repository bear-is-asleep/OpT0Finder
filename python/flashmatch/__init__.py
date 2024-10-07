import os
if not 'FMATCH_LIBDIR' in os.environ:
    raise ImportError('FMATCH_LIBDIR environment variable must be defined.')
from ROOT import gSystem
gSystem.Load(os.path.join(os.environ['FMATCH_LIBDIR'],'libflashmatch.so'))
from ROOT import flashmatch, phot, sim, geoalgo
#force loading C functions in dict by instantiating a class
c=flashmatch.PSet
c=flashmatch.load_pyutil
#from .demo import demo
from .toymc import ToyMC
from .rootinput import ROOTInput
from .utils import AnalysisManager