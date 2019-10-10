import json
import csv
import numpy as np
import pickle
import sqlite3

r = 3440 # radius of Earth in nautical miles

verbose = False
airports = {}

predef_ids  = [       'Rat',    'APPLE',    'LAKIE',    'ERORE',      'KSCR' ]
predef_lats = [  36.3894444,  40.556136,  40.829133,  40.954819,  35.7042745 ]
predef_lons = [ -76.9113889, -74.062253, -73.976792, -73.899233, -79.5042976 ]

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

print "Done building airport db."

def ll2xyz(lon, lat):
    return np.array( [np.cos(lat*np.pi/180) * np.cos(lon*np.pi/180),
                      np.cos(lat*np.pi/180) * np.sin(lon*np.pi/180),
                      np.sin(lat*np.pi/180)               ])

def llcoords(id):
    try:
        lat = airports[id]['lat']
        lon = airports[id]['lon']
    except KeyError:
        print "Can't find airport with ICAO ID '" + id + "'"
        raise
        return (0,0)
    return (lat, lon)

def ucoords(id):
    try:
        lat = airports[id]['lat']
        lon = airports[id]['lon']
    except KeyError:
        print "Can't find airport with ICAO ID '" + id + "'"
        raise
        return (0,0,0)
    return ll2xyz(lon, lat)

def dist(a,b):
    return r * np.abs(np.arccos(np.inner(a,b)))
    
def linedist(a,b,c):
    q = np.cross(a,b)
    n = q/np.sqrt(q.dot(q))
    return r * np.abs(np.arccos(np.inner(n,c))-np.pi/2)

log = []

'''
with open('log.csv') as logfile:
    logreader = csv.reader(logfile)
    header = logreader.next()
    deprt_i = header.index('Departure')
    destn_i = header.index('Destination')
    route_i = header.index('Route')
    durat_i = header.index('Duration')
    for row in logreader:
        flight = { 'deprt': row[deprt_i],
                   'destn': row[destn_i],
                   'route': row[route_i] if row[route_i] != '' else row[deprt_i]+'-'+row[destn_i],
                   'durat': float(row[durat_i]) }
        stops = flight['route'].split('-')
        flight['legs'] = [[stops[n], stops[n+1]] for n in range(len(stops)-1)]
        log.append(flight)
'''

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
    flight['legs'] = [[stops[n], stops[n+1]] for n in range(len(stops)-1)]
    log.append(flight)
    if verbose:
        print flight['route']
        print flight['legs']
db.close()

print "Done building list of flights."

nlon = 2**13
nlat = 2**12
lons = np.linspace(-180,180,nlon)
lats = np.linspace(-90,90,nlat)
lonv, latv = np.meshgrid(lons, lats, indexing='xy')
xyzs = ll2xyz(lonv,latv).transpose(1,2,0)
area_comp = np.cos(latv*np.pi/180)

def linemask(a,b,r=10):
    d = dist(a,b)
    return ( ( (linedist(a,b,xyzs)<r) &
               (dist(a,xyzs)<d)       &
               (dist(b,xyzs)<d)         )
             | (dist(a,xyzs)<r          )
             | (dist(b,xyzs)<r          ) )
             
def pointmask(a,r=15):
    return dist(a,xyzs)<r

bigmap = np.zeros((nlat,nlon))

print "Done building big blank map."

alllegs = []

for flight in log:
    littlemap = np.full((nlat,nlon),False,dtype='bool')
    for leg in flight['legs']:
        try:
            if leg[0] == leg[1]:
                littlemap |= pointmask(ucoords(leg[0]))
                alllegs.append([llcoords(leg[0])])
            else:
                littlemap |= linemask(ucoords(leg[0]),ucoords(leg[1]))
                alllegs.append([llcoords(leg[0]), llcoords(leg[1])])
        except KeyError:
            print "Problem with flight " + flight['route'] + "."
    area = (littlemap * area_comp).sum()
    if area > 0:
        bigmap += flight['durat'] * littlemap / area

print "Done plotting all the flights."

pickle.dump((bigmap, lats, lons, alllegs), open('log.p', 'wb'))

print "Done saving the map."




