Hi all,

We perform two corrections, one for temperature and one for flexure.  The temperature correction is -19.9 steps/degree F, so the focus position in steps would be reduced by ~20 for each Fahrenheit degree increase.

To find focus at a given position we typically take 7 images spaced by 25 focus steps at the 61".  The FWHM of the images are measured and plotted and a parabola fitted to the points.  The minimum gives best focus across the image. As an example, see 30imagefocus.png from the 1m using 30 images instead of 7.

The 61" flexure correction map is attached as a plot and a data file.  At the 61" focus does not change much with telescope position relative to the step range needed to find focus.  The map was generated from 78 points and smoothed. There is room for improvement but the corrections are small so it hasn't been a high priority.

focusmap.dat contains the correction map in an array of sorts, for easy lookup. The format is one long line containing:

{0.1 0.5 ...} {} {} ...

There are 91 groups of {}, one for each altitude from 0 to 90.  Within each {} is a list of 360 corrections for that altitude at every azimuth, 0 to 359. While this is generally overkill, especially at the 61", it allows for mapping details.  As you would expect, all the values in the last {}, which is alt 90, are the same since it's really one point.

Alex
