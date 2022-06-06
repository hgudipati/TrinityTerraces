# TrinityTerraces
Code for generating figures of my ESurf paper: A multi-proxy assessment of terrace formation in the lower Trinity River 
Valley, Texas. This file contains the code to load and plot plot the data for Figure 4 - 5, 7 - 12. 

## Input (data found in the supplementary): 

* S1 TrinityTerraces_characteristics.csv: include all terrace specific data including terrace id, classification, 
median and interquartile range of UTM northing, easting, elevation, and detrended elevation as well as 
terrace ids of adjacent lower and higher terraces. The minimum bounding box width and length are also 
included.
* S2: Terrace_5m_elevation_detrended_elevation.csv: includes all elevation and detrended elevations at 
the 5m DEM resolution and including their coordinates, corresponding terrace and terrace classification. 
* S3: Terrace_10m_valley_plane_elevation.csv: includes the corresponding elevations on terraces from 
the best-fit plane to the modern valley data.
* S4: Trinity_Terraces_Plane_Fitting_output.csv: includes the RMSE of a best-fit plane fitted to the high, 
intermediate and low Deweyville terraces as well as the RMSE of randomly classified terrace sets.
* S5: Paleochannel_summary.csv: includes all paleochannel specific data including paleochannel id, mean 
and standard deviation of width, length of paleochannel, corresponding terrace id and terrace 
classification, assigned representative grain size, number of channel bends preserved, and channel bend
appear to preserve bend cut-off.
* S6: Paleochannel_bends_on_terraces.csv: includes a closer look at the paleochannel bends preserved 
on each terrace and aggregates the total number of bends on each terraces, the maximum number of 
bends related to on generation of channel, and the maximum number of bends for one continuous 
channel segment. 

## Output:

* Figure 4: Median terrace elevation and (B) detrended elevation for the 52 terraces along the N-S trending valley, 
S7: PaleochannelWidthsMeasurements includes all width measurements for each paleo-channel used 
to find the mean and standard deviation of width as well as the length of the paleochannel in S5.
* Figure 5: Distributions of detrended elevations for terraces classified by Garvin (2008).  
* Figure 7. Root mean square error (RMSE) of a plane fit to elevation points of terraces previously classified as high Deweyville, intermediate Deweyville, and low Deweyville in the Trinity River valley compared to  distribution of RMSE from 150,000 randomly grouped terraces.
* Figure 8. Paleo-channel widths and paleo-discharge estimates. (A) Elevation transects for six paleo-channels (T1-T6). Transects are taken from locations indicated in (B)-(D) with mapped paleo-channels outlined in blue and terrace extents mapped outlined in grey. (E) Paleo-discharge estimates for the Trinity River are plotted as a function of their width. Each paleo-discharge was calculate using preserved channel width measurements and the discharge-width relationship from Wilkerson and Parker (2011) (Eq. 1). Error bars represent the first and third quartile of paleo-channel discharge estimates.
* Figure 9. Paleo-discharge estimates for the Trinity River plotted as a function of their associated detrended terrace elevations. 
* Figure 10. Mixing model fits to measured distributions of terrace elevation and estimated paleo-discharges. Distributions using (A) elevation and (B) paleo-channels support an interpretation of allogenic forcing for high terrace abandonment due to increasing discharge. Akaike Information Criterion (AIC) for (C) detrended elevation and (D) paleo-discharge mixing model.
* Figure 11. (A) terrace length, (B) paleo-channel length, and (C) paleo-channel width plotted along the median elevation of the associated terrace above the modern valley floor plane.
* Figure 12. Terrace properties used to assess the likelihood of meander bend-cutoff being the driver of terrace formation. (A) Differences in elevation between adjacent terrace surfaces. Also plotted as vertical lines and swaths are the mean values ± 1 standard deviation for elevation decreases expected from cutting off a single meander bend for paleo-channels of the low, intermediate, and high Deweyville bounding surfaces. (B) Maximum number of paleo-meander bends preserved in a channel segment on each terrace. Most terraces have between 0-1 channel bends preserved for one generation of channel.
