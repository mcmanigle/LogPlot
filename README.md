# LogPlot

This is a simple pair of Python scripts to take a pilot's logbook and format it as a basic heatmap.  Written in response to the various fill-in-the-state maps on the Blue Board.

The first script, `logmap.py`, reads a CSV pilot log file (`log.csv`) and outputs a pickle archive, `log.p`, which contains the basic map data.  The second script, `logplot.py`, uses `log.p` and outputs PNG graphics displaying the heatmap.

This isn't in any way production-ready software, so some tweaking will undoubtedly be required.  It requires the [basemap] (http://matplotlib.org/basemap/) package, as well as numpy and scipy.  Out of the box, the pilot log file should be CSV formatted, with headers in the first row, and at least the following columns:
* `Departure` - ICAO identifier of departure airport.
* `Destination` - ICAO identifier of destination airport.
* `Duration` - Time en route, formatted as a float.
* `Route` - Series of ICAO airport identifiers, separated by dashes (`-`), describing route (may be blank).

## Example Images

Included in this package is an example flight log and resulting images:

![First heatmap example]
(1.png)

![Second heatmap example]
(2.png)