import pickle
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import scipy.ndimage
import scipy.misc
import matplotlib as mpl
from PIL import Image
mpl.use('agg')

bigmap, lats, lons, alllegs = pickle.load( open('log.p', 'rb') )
lonv, latv = np.meshgrid(lons, lats, indexing='xy')

print('Done loading map data.')

bigmap = bigmap.toarray()
bigmap /= bigmap[np.nonzero(bigmap)].min()
bigmap = np.power(bigmap, .3)
maxz = np.max(bigmap)


r  = 6378100 # radius of earth in meters

def blockshaped(arr, nrows, ncols):
    """
    Return an array of shape (n, nrows, ncols) where
    n * nrows * ncols = arr.size

    If arr is a 2D array, the returned array should look like n subblocks with
    each subblock preserving the "physical" layout of arr.
    """
    h, w = arr.shape
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def paramline(p1, p2, ps=25):
    xs = np.linspace(p1[0], p2[0], ps)
    ys = np.linspace(p1[1], p2[1], ps)
    return (xs, ys)
    
def paramcurve(p1, p2, offsetlim=0.2, ps=25):
    xs, ys = paramline(p1, p2, ps)
    t = np.sin(np.linspace(0, np.pi, ps))
    d = np.linalg.norm(np.array(p1)-np.array(p2))
    if offsetlim > d * 0.2:
        offsetlim = d * 0.2
    offset = np.random.uniform(-offsetlim, offsetlim)
    xo = np.array(offset * (p2[1]-p1[1])/d)
    yo = np.array(offset * (p2[0]-p1[0])/d)
    return (xs + xo*t, ys + yo*t)

def localcircle(p, r, theta=0, ps=25):
    t = np.linspace(0, 2*np.pi, ps)
    c = (r * np.cos(theta), r*np.sin(theta))
    xs = p[0] + c[0] + r * np.cos(t)
    ys = p[1] + c[1] + r * np.sin(t)
    return (xs, ys)

def localpath_old(p, r, theta=0, ps=25):
    t = np.linspace(0, 2*np.pi, ps)
    xs = r * np.cos(t)
    ys = 0.6 * r * np.sin(2*t)
    return (xs * np.cos(theta) - ys * np.sin(theta) + p[0],
            xs * np.sin(theta) + ys * np.cos(theta) + p[1] )

def localpath(p, ps=25):
    radius = np.random.uniform(0.1, 0.2)
    th_0 = np.random.uniform(0, 2.*np.pi)
    th_d = np.random.uniform(np.pi, 7.*np.pi/4.)
    r_sin_deltpct = np.random.uniform(5., 10.)
    r_sin_deltn = np.random.choice([4,5])
    r_lin_deltpct = np.random.uniform(-30., 30.)

    tr  = np.linspace(0, np.pi, ps)
    r   = radius * ( np.sin(tr) +
                     r_sin_deltpct * np.sin(tr * (1 + r_sin_deltn * 2)) / 100. +
                     r_lin_deltpct * np.sin(tr) * tr / ( 100. * np.pi )
                   )
    th  = np.linspace(th_0, th_0 + th_d, ps)
    return ( p[0] + r * np.cos(th), p[1] + r * np.sin(th) )

lonn = 16
latn = 8
bmgrid = blockshaped(bigmap, bigmap.shape[0]//latn, bigmap.shape[1]//lonn)
lonvgr = blockshaped(lonv,   bigmap.shape[0]//latn, bigmap.shape[1]//lonn)
latvgr = blockshaped(latv,   bigmap.shape[0]//latn, bigmap.shape[1]//lonn)

bmgridtf = np.array([x.sum()>0 for x in bmgrid])

gridbits = []

adj_x = [ -1, -1, -1,  0,  0,  1,  1,  1 ]
adj_y = [ -1,  0,  1, -1,  1, -1,  0,  1 ]

def adding_to_gridbit(bit, gridtf, dims):
    refbit = bit.copy()
    adding = False
    for i in refbit:
        for delta in zip(adj_x, adj_y):
            j = np.array(delta) + np.unravel_index(i, dims)
            if gridtf[np.ravel_multi_index(j, dims, 'wrap')]:
                gridtf[np.ravel_multi_index(j, dims, 'wrap')] = False
                bit.add(np.ravel_multi_index(j, dims, 'wrap'))
                adding = True
    return adding

for i in range(bmgridtf.size):
    if not bmgridtf[i]: continue
    gridbit = set([i])
    bmgridtf[i] = False
    adding = True
    while adding:
        adding = adding_to_gridbit(gridbit, bmgridtf, (latn, lonn))
    gridbits.append(gridbit)

submaps = []
border = 0.2

for b in gridbits:
    lonmin, lonmax, latmin, latmax = (np.inf, -np.inf, np.inf, -np.inf)
    xs, ys, zs = (np.array([]), np.array([]), np.array([]))
    for i in b:
        lonmask = np.ma.masked_array(lonvgr[i], mask = bmgrid[i] < 0.5)
        latmask = np.ma.masked_array(latvgr[i], mask = bmgrid[i] < 0.5)
        lonmin = np.min([lonmin, np.ma.min(lonmask)])
        lonmax = np.max([lonmax, np.ma.max(lonmask)])
        latmin = np.min([latmin, np.ma.min(latmask)])
        latmax = np.max([latmax, np.ma.max(latmask)])
    
    lonspan = lonmax - lonmin
    latspan = latmax - latmin
    lonmax += border * lonspan
    lonmin -= border * lonspan
    latmax += border * latspan
    latmin -= border * latspan
        
    submap = { 'lonmin': lonmin,
               'lonmax': lonmax,
               'latmin': latmin,
               'latmax': latmax,
               'lonis':  np.searchsorted(lons, [lonmin, lonmax]),
               'latis':  np.searchsorted(lats, [latmin, latmax])
             }
    submaps.append(submap)

print('Done generating submap data.')

from matplotlib import cm

cmap_resolution = 100
cmap = cm.get_cmap('jet', cmap_resolution)
cmap_vals = cmap(np.arange(cmap_resolution)) #extract those values as an array 
fade_size = 12
first_stage_v = 0.6
first_stage_z = 3*cmap_resolution//10
second_stage_v = 0.9
second_stage_z = 6*cmap_resolution//10

cmap_vals[0][3] = 0.0 #change the first value 
for i in range(1, fade_size):
    cmap_vals[i][3] = first_stage_v * ((i-1.0)/fade_size)**1.5
for i in range(fade_size, first_stage_z):
    cmap_vals[i][3] = first_stage_v
for i in range(first_stage_z, second_stage_z):
    cmap_vals[i][3] = first_stage_v + (second_stage_v - first_stage_v) * (i-first_stage_z)/(second_stage_z - first_stage_z)
for i in range(second_stage_z, cmap_resolution):
    cmap_vals[i][3] = second_stage_v

opcmap = mpl.colors.LinearSegmentedColormap.from_list('opcmap', cmap_vals)


import scipy.ndimage
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from cartopy.io.img_tiles import Stamen

def figsize(w, h, figscale = 25):
    if w > h:
        return (figscale, h*figscale/w)
    else:
        return (w*figscale/h, figscale)

projection = ccrs.PlateCarree()
dpi = 300
zoom = 5
blur = 24
quick = False

for i, submap in enumerate(submaps):
    
    lonmid = (submap['lonmax']+submap['lonmin'])/2
    latmid = (submap['latmax']+submap['latmin'])/2
    width  = r * (submap['lonmax']-submap['lonmin']) * np.pi / 180 * np.cos( latmid * np.pi / 180 )
    height = r * (submap['latmax']-submap['latmin']) * np.pi / 180
    
    if quick:
        lons = lonv[(submap['latis'][0]):(submap['latis'][1]),
                    (submap['lonis'][0]):(submap['lonis'][1])]
        lats = latv[(submap['latis'][0]):(submap['latis'][1]),
                    (submap['lonis'][0]):(submap['lonis'][1])]
        zs   = bigmap[(submap['latis'][0]):(submap['latis'][1]),
                      (submap['lonis'][0]):(submap['lonis'][1])]
    
    else:
        lons = scipy.ndimage.interpolation.zoom(
                                  lonv[(submap['latis'][0]):(submap['latis'][1]),
                                       (submap['lonis'][0]):(submap['lonis'][1])],
                                  zoom, order=1)
        lats = scipy.ndimage.interpolation.zoom(
                                  latv[(submap['latis'][0]):(submap['latis'][1]),
                                       (submap['lonis'][0]):(submap['lonis'][1])],
                                  zoom, order=1)
        zs   = scipy.ndimage.interpolation.zoom(
                                  bigmap[(submap['latis'][0]):(submap['latis'][1]),
                                         (submap['lonis'][0]):(submap['lonis'][1])],
                                  zoom, order=3)
        zs = gaussian_filter(zs, sigma=blur)
    
    print('Done zooming and blurring map data for submap ' + str(i+1) + '.')
    

    f = plt.figure(figsize=figsize(width, height), dpi=dpi, tight_layout=True)
    logmap = plt.axes(projection=projection)

    tiler = Stamen('terrain-background')
    logmap.add_image(tiler, 10, interpolation='spline36')
    
    state_boundaries = cfeature.NaturalEarthFeature('cultural', 'admin_1_states_provinces', '10m')
    logmap.add_feature(state_boundaries, facecolor=(0,0,0,0), edgecolor=(0.2,0.2,0.2,0.5))

    print('Done plotting base map for submap ' + str(i+1) + '.')

    for leg in alllegs:
        if len(leg) == 1:
            p = (leg[0][1], leg[0][0])
            xs, ys = localpath(p)
        elif len(leg) == 2:
            p1 = (leg[0][1], leg[0][0])
            p2 = (leg[1][1], leg[1][0])
            xs, ys = paramcurve(p1, p2)
        else:
            continue
        logmap.plot(xs, ys, linewidth=1., color='white')
    
    print('Done plotting paths for submap ' + str(i+1) + '.')

    logmap.pcolormesh(lons,lats,zs,cmap=opcmap,vmax=maxz, zorder=100)

    print('Done plotting color on submap ' + str(i+1) + '.')

    logmap.set_extent([submap['lonmin'], submap['lonmax'],
                       submap['latmin'], submap['latmax']], crs=projection)
    plt.savefig(str(i+1)+'.png', dpi='figure')
    
    print('Saved submap ' + str(i+1) + '.')
