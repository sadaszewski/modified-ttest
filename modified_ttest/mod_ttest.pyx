from libc.stdint cimport int32_t, \
    uintptr_t
from libc.string cimport memcpy
from libc.math cimport sqrt
import cython

import numpy as np
cimport numpy as np
np.import_array()

from libc.stdio cimport printf


cdef extern:
    struct DIMS_struct:
        int n
        int p
        int nclass
        
    double distance_max(double *xpos, double *ypos, int n)
    
    void set_bounds(DIMS_struct *dims, double maxdist, int do_half, double *upper_bounds)
    
    void mod_ttest(double *x, double *y, DIMS_struct *dims,
                   double *xpos, double *ypos, double *upper_bounds,
                   double *cor, double *card, double *imoran, double *stats)
    

def modified_ttest(x, y, coords, nclass=13):
    cdef np.double_t[::1] x_buf
    cdef np.double_t[::1] y_buf
    cdef np.double_t[::1] xpos_buf
    cdef np.double_t[::1] ypos_buf
    cdef np.int32_t[::1] dims_buf
    cdef np.double_t[:, ::1] corr_buf
    cdef np.double_t[::1] upper_bounds_buf
    cdef np.double_t[::1] card_buf
    cdef np.double_t[::1] imoran_buf
    cdef np.double_t[::1] stats_buf
    
    if len(x) != len(y):
        raise ValueError('x and y must have the same length')
    x = np.array(x, dtype=np.double)
    y = np.array(y, dtype=np.double)
    coords = np.array(coords, dtype=np.double)
    n = len(x)
    ndist = n * (n - 1) / 2
    if nclass is None:
        nclass = int(1.5 + 3.3 * np.log10(ndist))
    xpos = coords[:, 0].copy()
    ypos = coords[:, 1].copy()
    p = coords.shape[1]
    if p != 2:
        raise ValueError('Only implemented for 2D coords')
    dims = np.array([n, p, nclass], dtype=np.int32)
    corr = np.corrcoef(x, y)
    
    upper_bounds = np.zeros(nclass, dtype=np.double)
    card = np.zeros(nclass, dtype=np.double)
    imoran = np.zeros(nclass * p, dtype=np.double)
    stats = np.zeros(3, dtype=np.double)
    
    x_buf = x
    y_buf = y
    xpos_buf = xpos
    ypos_buf = ypos
    dims_buf = dims
    corr_buf = corr
    upper_bounds_buf = upper_bounds
    card_buf = card
    imoran_buf = imoran
    stats_buf = stats
    
    maxdist = distance_max(&xpos_buf[0], &ypos_buf[0], n)
    do_half = 0
    set_bounds(<DIMS_struct*> &dims_buf[0], maxdist, do_half, &upper_bounds_buf[0])

    mod_ttest(
        &x_buf[0],
        &y_buf[0],
        <DIMS_struct*> &dims_buf[0],
        &xpos_buf[0],
        &ypos_buf[0],
        &upper_bounds_buf[0],
        &corr_buf[0, 1],
        &card_buf[0],
        &imoran_buf[0],
        &stats_buf[0])
    
    return stats
