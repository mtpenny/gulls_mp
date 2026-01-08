import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys


from matplotlib.collections import LineCollection


def colored_line(x, y, c, ax, **lc_kwargs):
    """
    Plot a line with a color specified along the line by a third value.

    It does this by creating a collection of line segments. Each line segment is
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should be the same size as x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """
    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Default the capstyle to butt so that the line segments smoothly line up
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)

    # Compute the midpoints of the line segments. Include the first and last points
    # twice so we don't need any special syntax later to handle them.
    x = np.asarray(x)
    y = np.asarray(y)
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    # Determine the start, middle, and end coordinate pair of each line segment.
    # Use the reshape to add an extra dimension so each pair of points is in its
    # own list. Then concatenate them to create:
    # [
    #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
    #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
    #   ...
    # ]
    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    lc = LineCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment

    return ax.add_collection(lc)



if len(sys.argv)==1:
    print(f"Usage: python {sys.argv[0]} <lightcurve>")
    exit()

lightcurve = sys.argv[1]

data = pd.read_csv(lightcurve,sep=r'\s+',comment='#')
outfile = lightcurve[:lightcurve.rfind('_')] +'.out'
idx = lightcurve[lightcurve.rfind('_')+1:lightcurve.find('.')]
print(lightcurve,outfile,idx)
out = pd.read_csv(outfile,sep='\s+')
outdata = out[out['EventID']==int(idx)].squeeze()
print(outdata)

nlens = pd.Series(list(data.columns)).str.contains('lens').sum()//2
print(f"nlens = {nlens}")

#header = pd.read_csv(sys.argv[1],sep='\s+',header=None,comment=None,engine='python',nrows=50,index_col=False)
#fsm = header[header.iloc[:,0]=='#fs:'].squeeze(axis=0)[1:].astype(float)
#event = header[header.iloc[:,0]=='#Event:'].squeeze(axis=0)[1:].astype(float)
#planet = header[header.iloc[:,0]=='#Planet:'].squeeze(axis=0)[1:].astype(float)
#source = header[header.iloc[:,0]=='#Obssrcmag:'].squeeze(axis=0)[1:].astype(float)
#lens = header[header.iloc[:,0]=='#Obslensmag:'].squeeze(axis=0)[1:].astype(float)

ls = ['-','--','-.']

plt.figure()

for i in range(nlens):
    plt.plot(data[f"lens{i}_x"],data[f"lens{i}_y"],label=f'{i}')

plt.gca().set_aspect('equal')
plt.legend()
#plt.colorbar(label='Time [days]')
plt.xlabel(r'$x$ [$r_{\rm E}$]')
plt.xlabel(r'$y$ [$r_{\rm E}$]')
plt.show()
