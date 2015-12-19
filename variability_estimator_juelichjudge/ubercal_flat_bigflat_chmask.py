import numpy
import calib


night_zero = 54900.
cellperflat = 1
nfterm = (64/cellperflat)**2

def nindx_mjd_fun(mjd):
    return numpy.uint32(numpy.floor(mjd)-night_zero)

flat_seasons = [nindx_mjd_fun(f) for f in [54900, 55296, 55327, 55662, 56110]]


def unique(obj):
    """Gives an array indexing the unique elements in the sorted list obj.

    Returns the list of indices of the last elements in each group
    of equal-comparing items in obj
    """
    nobj = len(obj)
    if nobj == 0:
        return numpy.zeros(0, dtype='i8')
    if nobj == 1:
        return numpy.zeros(1, dtype='i8')
    out = numpy.zeros(nobj, dtype=numpy.bool)
    out[0:nobj-1] = (obj[0:nobj-1] != obj[1:nobj])
    out[nobj-1] = True
    return numpy.sort(numpy.flatnonzero(out))

class subslices:
    "Iterator for looping over subsets of an array"
    def __init__(self, data):
        self.uind = unique(data)
        self.ind = 0
    def __iter__(self):
        return self
    def __len__(self):
        return len(self.uind)
    def next(self):
        if self.ind == len(self.uind):
            raise StopIteration
        if self.ind == 0:
            first = 0
        else:
            first = self.uind[self.ind-1]+1
        last = self.uind[self.ind]+1
        self.ind += 1
        return first, last

def apply_flat_solution(obj, soldict, include_resid=False):
    s = numpy.argsort(obj['filterid'])
    zp = numpy.zeros(len(obj), dtype='f4')
    zp[:] = numpy.nan
    findx = query_to_findx(obj)
    for f,l in subslices(obj['filterid'][s]):
        ind = s[f:l]
        filter = obj[ind[0]]['filterid'][0]
        sol0 = soldict[filter]
        tzp, tflat = sol0[0:2]
        if len(sol0) == 3:
            chmask = sol0[2]
        mind = numpy.searchsorted(tzp['mjd_obs'], obj['mjd_obs'][ind])
        m = ((mind >= 0) & (mind < len(tzp)) & (findx[ind] >= 0) &
             (findx[ind] < len(tflat)))
        m[m] &= (tzp['mjd_obs'][mind[m]] == obj['mjd_obs'][ind[m]])
        zp[ind[m]] = numpy.array(tzp[mind[m]]['zp'] + tflat[findx[ind[m]]],
                                 dtype='f4')
        if include_resid:
            zp[ind[m]] += numpy.array(tzp[mind[m]]['resid'], dtype='f4')
        if len(sol0) == 3:
            badchip = chmask['chmask'][mind[m], obj['chip_id'][ind[m]]] != 0
            zp[ind[m][badchip]] = numpy.nan
    return zp

def query_to_findx(query):
    cellid = calib.ps_xy2cell(query['x_psf'], query['y_psf'],
                              chip_id=query['chip_id'])
    xim, yim = calib.cell2flat_ij(query['chip_id'], cellid)
    xim /= cellperflat ; yim /= cellperflat
    out = (yim*(64/cellperflat)+xim).astype('i4')
    if flat_seasons is not None:
        nperflat = (64/cellperflat)**2
        ind = numpy.searchsorted(flat_seasons, query['mjd_obs']-54900)
#        if numpy.any(ind <= 0):
#            pdb.set_trace()
        out = out + nperflat * (ind - 1)
    return out

def read_flat_solution(filename, chmask=False):
    import pyfits
    filters = 'grizy'
    sol = { }
    increment = 2 + (1 if chmask else 0)
    for i,f in enumerate(filters):
        zp = pyfits.getdata(filename, increment*i+1)
        flat = pyfits.getdata(filename, increment*i+2)
        sol0 = zp, flat
        if chmask:
            chmask0 = pyfits.getdata(filename, increment*i+3)
            sol0 = sol0 + (chmask0,)
        sol[f] = sol0
    return sol
