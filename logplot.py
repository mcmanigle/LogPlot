import pickle
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import scipy.ndimage
import matplotlib as mpl

bigmap, lats, lons, alllegs = pickle.load( open('log.p', 'rb') )
lonv, latv = np.meshgrid(lons, lats, indexing='xy')

print 'Done loading map data.'

bigmap /= bigmap[np.nonzero(bigmap)].min()
bigmap = np.power(bigmap, .41)
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

def paramline(p1, p2, ps=100):
    xs = np.linspace(p1[0], p2[0], ps)
    ys = np.linspace(p1[1], p2[1], ps)
    return (xs, ys)
    
def paramcurve(p1, p2, offsetpct, ps=50):
    xs, ys = paramline(p1, p2, ps)
    t = np.sin(np.linspace(0, np.pi, ps))
    d = np.linalg.norm(np.array(p1)-np.array(p2))
    offset = 0.01 * offsetpct * d
    xo = np.array(offset * (p2[1]-p1[1])/d)
    yo = np.array(offset * (p2[0]-p1[0])/d)
    return (xs + xo*t, ys + yo*t)

lonn = 16
latn = 8
bmgrid = blockshaped(bigmap, bigmap.shape[0]/latn, bigmap.shape[1]/lonn)
lonvgr = blockshaped(lonv,   bigmap.shape[0]/latn, bigmap.shape[1]/lonn)
latvgr = blockshaped(latv,   bigmap.shape[0]/latn, bigmap.shape[1]/lonn)

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

print 'Done generating submap data.'

from matplotlib import cm

cmap_resolution = 50
cmap = cm.get_cmap('jet', cmap_resolution) #generate a jet map with 10 values 
cmap_vals = cmap(np.arange(cmap_resolution)) #extract those values as an array 
cmap_vals[0][3] = 0.0 #change the first value 
fade_size = 15
for i in range(1, fade_size+1):
    cmap_vals[i][3] = ((i-1.0)/fade_size)**1.5

newcmap = mpl.colors.LinearSegmentedColormap.from_list('newcmap', cmap_vals)

import scipy.ndimage
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

figsize = (15,15)
zoom = 20
blur = 30
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

    plt.figure(figsize=figsize)
    map = Basemap(width=width, height=height,lon_0=lonmid,lat_0=latmid,
                resolution='i',projection='cass')

    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.drawstates(linewidth=0.25)
    #map.fillcontinents(color='coral',lake_color='aqua')
    #map.drawmapboundary(fill_color='aqua')
    #map.drawmeridians(np.arange(0,360,30))
    #map.drawparallels(np.arange(-90,90,30))
    pc = map.pcolormesh(lons,lats,zs,latlon=True,vmax=maxz,cmap=newcmap)
    
    for leg in alllegs:
        if len(leg) <= 1:
            continue
        offset = np.random.uniform(-5,5)
        p1 = map(leg[0][1], leg[0][0])
        p2 = map(leg[1][1], leg[1][0])
        xs, ys = paramcurve(p1, p2, offset)
        map.plot(xs, ys, linewidth=1.8, color=(1,1,1,0.2))

    plt.savefig(str(i+1)+'.png')
    plt.close()

