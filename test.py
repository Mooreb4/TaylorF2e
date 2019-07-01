import ctypes, numpy, pycbc.conversions
from pycbc.types import FrequencySeries
from ctypes import c_double, c_void_p

lib = ctypes.cdll.LoadLibrary("./tf2e.so")

fend = 1000.0
df = 1.0 / 256
flen = int(fend / df)
flow = 30.0
lamc = 0.
lc = 0.
dist = 1e18
inc = 0.1
ecc = 0.2
m1 = m2 = 1.4

mchirp = float(pycbc.conversions.mchirp_from_mass1_mass2(m1, m2))
eta = float(pycbc.conversions.eta_from_mass1_mass2(m1, m2))

hp = numpy.zeros(flen, dtype=numpy.complex128)
hc = hp.copy()

f = lib.generate
f.argtypes = [c_void_p, c_void_p, c_double, c_double, c_double,
              c_double, c_double, c_double, c_double, c_double, c_double]
_ = f(hp.ctypes.data, hc.ctypes.data,
             mchirp, eta, inc, ecc, lamc, lc, dist, fend, df)

import pylab
# it appears that the plus / cross data is time inverted
hp = FrequencySeries(hp.conj(), delta_f=df, epoch=-int(1.0 / df))
hc = FrequencySeries(hc.conj(), delta_f=df, epoch=-int(1.0 / df))

kmin = int(flow / df)
hp[:kmin].clear()
hc[:kmin].clear()

t = hp.to_timeseries()

pylab.plot(t.sample_times, t)
#pylab.xscale('log')
pylab.show()
