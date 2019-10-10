import json
import csv
import numpy as np
from scipy import sparse
import pickle
import sqlite3
import multiprocessing
from Queue import Empty
from functools import partial

verbose = False

# r = 3440 # radius of Earth in nautical miles

def ll2xyz(lon, lat):
    return np.array( [np.cos(lat*np.pi/180) * np.cos(lon*np.pi/180),
                      np.cos(lat*np.pi/180) * np.sin(lon*np.pi/180),
                      np.sin(lat*np.pi/180)               ])

def dist(a,b):
    return 3440 * np.abs(np.arccos(np.inner(a,b)))
    
def linedist(a,b,c):
    q = np.cross(a,b)
    n = q/np.sqrt(q.dot(q))
    return 3440 * np.abs(np.arccos(np.inner(n,c))-np.pi/2)

def calculate_map_coordinates():
    nlon = 2**14
    nlat = 2**13

    lons = np.linspace(-180,180,nlon)
    lats = np.linspace(-90,90,nlat)
    lonv, latv = np.meshgrid(lons, lats, indexing='xy')
    xyzs = ll2xyz(lonv,latv).transpose(1,2,0)
    print("Done calculating map coordinates.")
    return xyzs, lats, lons
    
def build_airport_db():
    airports = {}

    predef_ids  = [      'NC18',       'Rat',    'APPLE',    'LAKIE',    'ERORE',      'KSCR' ]
    predef_lats = [  36.3894444,  36.3894444,  40.556136,  40.829133,  40.954819,  35.7042745 ]
    predef_lons = [ -76.9113889, -76.9113889, -74.062253, -73.976792, -73.899233, -79.5042976 ]

    for id, lat, lon in zip(predef_ids, predef_lats, predef_lons):
        ap = { 'lat': lat, 'lon': lon }
        airports[id] = ap

    with open('airports.csv') as airportdb:
        apreader = csv.reader(airportdb)
        header = apreader.next()
        idnt_i = header.index('ident')
        iata_i = header.index('iata_code')
        locl_i = header.index('local_code')
        lat_i  = header.index('latitude_deg')
        lon_i  = header.index('longitude_deg')
        for row in apreader:
            ap = { 'lat': float(row[lat_i]),
                   'lon': float(row[lon_i])  }
            airports[row[idnt_i]] = ap
            if row[iata_i] not in airports:
                airports[row[iata_i]] = ap
            if row[locl_i] not in airports:
                airports[row[locl_i]] = ap

    print("Done building airport db.")
    return airports

def build_flight_log():
    log = []

    db = sqlite3.connect('/Users/mcmanigle/Dropbox/Apps/PilotPro/Logbook.pilotpro')

    route_cfid = db.execute('SELECT customFieldId from customFields WHERE name="Route"').fetchone()[0]

    flights = db.execute('SELECT l.departure, l.destination, c.value, l.duration ' +
                         'FROM logbookEntries AS l LEFT JOIN customValues AS c ' +
                         'ON l.logbookEntryId = c.logbookEntryId AND c.customFieldId = "' + route_cfid + '"')
                     
    for row in flights:
        flight = { 'deprt': row[0],
                   'destn': row[1],
                   'route': row[2] if row[2] else row[0]+'-'+row[1],
                   'durat': float(row[3]) }
        stops = flight['route'].split('-')
        if( stops[0]  != row[0] ): stops = [row[0]] + stops
        if( stops[-1] != row[1] ): stops += [row[1]]
        flight['legs'] = [[stops[n], stops[n+1]] for n in range(len(stops)-1)]
        log.append(flight)
        if verbose:
            print(flight['route'])
            print(flight['legs'])
    db.close()
    
    print("Done building list of flights.  " + str(len(log)) + " total.")
    return log

def maps_legs_from_flights(process_num, log, airports, xyzs, bigmapq, alllegsq):
    nlat, nlon, _ = xyzs.shape
    bigmap = sparse.bsr_matrix((nlat,nlon), dtype='f4')

    max_degree_length = 60 # approx max lat/lon one degree distance in nautical miles
    lon_gridpt_length = max_degree_length * 360 / nlon  # approximate size of lat/lon
    lat_gridpt_length = max_degree_length * 180 / nlat  #   grid at equator
    
    lat_to_gridpt = lambda lat: (lat +  90) * nlat / 180
    lon_to_gridpt = lambda lon: (lon + 180) * nlon / 360
    
    lat_area_correct = lambda g: np.cos(g * np.pi / nlat - np.pi/2)
        
    def llcoords(id):
        try:
            lat = airports[id]['lat']
            lon = airports[id]['lon']
        except KeyError:
            print("  * Can't find airport with ICAO ID '" + id + "'")
            raise
            return (0,0)
        return (lat, lon)

    def ucoords(id):
        try:
            lat = airports[id]['lat']
            lon = airports[id]['lon']
        except KeyError:
            print("  * Can't find airport with ICAO ID '" + id + "'")
            raise
            return (0,0,0)
        return ll2xyz(lon, lat)
    
    def linemask(map, a,b,r=10):
        ua = ll2xyz(a[1],a[0])
        ub = ll2xyz(b[1],b[0])
        d = dist(ua,ub)
        
        mingridlon = np.min(lon_to_gridpt(np.array([a[1], b[1]]))) - r/lon_gridpt_length
        maxgridlon = np.max(lon_to_gridpt(np.array([a[1], b[1]]))) + r/lon_gridpt_length
        mingridlat = np.min(lat_to_gridpt(np.array([a[0], b[0]]))) - r/lat_gridpt_length
        maxgridlat = np.max(lat_to_gridpt(np.array([a[0], b[0]]))) + r/lat_gridpt_length
        
        for glat in range(int(mingridlat), int(maxgridlat)):
            for glon in range(int(mingridlon), int(maxgridlon)):
                map[glat, glon] += ( ( (linedist(ua,ub,xyzs[glat, glon])<r) &
                                        (dist(ua,xyzs[glat, glon])<d)       &
                                        (dist(ub,xyzs[glat, glon])<d)         )
                                      | (dist(ua,xyzs[glat, glon])<r          )
                                      | (dist(ub,xyzs[glat, glon])<r          ) )
             
    def pointmask(map, a,r=15):
        mingridlon = lon_to_gridpt(a[1]) - r/lon_gridpt_length
        maxgridlon = lon_to_gridpt(a[1]) + r/lon_gridpt_length
        mingridlat = lat_to_gridpt(a[0]) - r/lat_gridpt_length
        maxgridlat = lat_to_gridpt(a[0]) + r/lat_gridpt_length

        for glat in range(int(mingridlat), int(maxgridlat)):
            for glon in range(int(mingridlon), int(maxgridlon)):
                map[glat, glon] += (dist(ll2xyz(a[1],a[0]),xyzs[glat, glon])<r)

    counter = 0
    for flight in log:
        littlemap = sparse.dok_matrix((nlat,nlon),dtype='bool')
        for leg in flight['legs']:
            try:
                if leg[0] == leg[1]:
                    pointmask(littlemap, llcoords(leg[0]))
                    alllegsq.put([llcoords(leg[0])])
                    alllegsq.cancel_join_thread()
                else:
                    linemask(littlemap, llcoords(leg[0]),llcoords(leg[1]))
                    alllegsq.put([llcoords(leg[0]), llcoords(leg[1])])
                    alllegsq.cancel_join_thread()
            except KeyError:
                print(" - Problem with flight " + flight['route'] + ".")
        area = 0
        for (lat, lon) in littlemap.keys():
            area += lat_area_correct(lat) * littlemap[lat, lon]
        if area > 0:
            bigmap += flight['durat'] * littlemap / area
        counter += 1
        if counter % 5 == 0: print(' - Process ' + str(process_num) + ' finished ' + str(counter) + ' flights.')
    
    bigmapq.put(bigmap)
    print('Process ' + str(process_num) + ' complete after ' + str(counter) + ' flights.')

if __name__ == '__main__':
    airports = build_airport_db()
    
    xyzs, lats, lons = calculate_map_coordinates()

    log = build_flight_log()

    num_jobs = multiprocessing.cpu_count() - 3
    bigmapq = multiprocessing.Queue()
    alllegsq = multiprocessing.Queue()
    
    print("Starting " + str(num_jobs) + " workers, for about " + str(len(log)/num_jobs) + " flights per worker...")
    
    for i in range(num_jobs):
        j = multiprocessing.Process( target=maps_legs_from_flights, 
                                     args=(i+1, log[i::num_jobs], airports, xyzs, bigmapq, alllegsq) )
        j.start()

    nlat, nlon, _ = xyzs.shape
    bigmap = sparse.bsr_matrix((nlat,nlon), dtype='f4')
    alllegs = []
    
    for i in range(num_jobs):
        bigmap += bigmapq.get()
    
    while True:
        try:           alllegs.append(alllegsq.get_nowait())
        except Empty:  break

    print("Done plotting all the flights.")

    pickle.dump((bigmap, lats, lons, alllegs), open('log.p', 'wb'))

    print("Done saving the map.")

