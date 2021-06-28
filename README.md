# LogPlot

This is a simple pair of Python scripts to take a pilot's logbook and format it as a basic heatmap.  Written in response to the various fill-in-the-state maps on the Blue Board.

The first script, `logmap.py`, reads a CSV pilot log file (`log.csv`) and outputs a pickle archive, `log.p`, which contains the basic map data.  The second script, `logplot.py`, uses `log.p` and outputs PNG graphics displaying the heatmap.

This isn't in any way production-ready software, so some tweaking will undoubtedly be required.  It requires the [cartopy] (https://scitools.org.uk/cartopy/) package, as well as [numpy] (https://numpy.org) and [scipy] (https://www.scipy.org), which can all be installed via [conda] (https://docs.conda.io/en/latest/).  Out of the box, the pilot log file is expected to be a [sqlite] (https://www.sqlite.org/) database like the one saved as a backup of [Pilot Pro] (https://pilotpro.com/) logbook (version 3, which is what I use). A little bit of tweaking would make it usable for any log format that includes:
* `Departure` - ICAO identifier of departure airport.
* `Destination` - ICAO identifier of destination airport.
* `Duration` - Time en route, formatted as a float.
* `Route` - Series of ICAO airport identifiers, separated by dashes (`-`), describing route (may be blank).

## Example Images

![First heatmap example]
(1.png)
