# Image Analysis Scripts

MATLAB scripts for EVL image analysis. Require MATLAB R2020b or later.

---

## EVL_surface_projection.m

Generates apical surface-projected images from a confocal z-stack TIFF file.

**What it does:** Detects the EVL apical surface using a reference fluorescence channel (e.g. membrane marker), then max-projects a target channel (e.g. Myosin II, Kibra) within a defined z-window around that surface. This isolates apical signals from basal/lateral contributions.

**Requires:** Image Processing Toolbox

**Key parameters to set before running:**

| Parameter | Description |
|---|---|
| `File_name` | Input TIFF file (multi-channel z-stack) |
| `Surface_ch` | Channel index used to detect apical surface |
| `Target_ch` | Channel index to project |
| `offset` | Z-offset from detected surface (slices) |
| `DeltaZ` | Half-width of projection window (slices) |
| `Th_value` | Intensity threshold for surface detection |
| `G_blur` | Gaussian filter size for smoothing (0 = skip) |
| `t_size, z_size, c_size` | TIFF dimensions (time, z, channel) |

**Outputs:**
- `Target_surface_max_*.tif` — projected target channel
- `Reference_surface_max_*.tif` — projected reference channel
- `Hight_map.tif` — detected surface z-position map

---

## Actin_kibra_xcorr.m

Computes spatial cross-correlation between Actin and Kibra fluorescence signals within each segmented EVL cell. Used to quantify the spatial association between Kibra condensates and actin-based apical projections.

**What it does:** For each cell in a Cellpose segmentation mask, erodes the mask to exclude junctional signals, computes 2D normalized cross-correlation between Actin and Kibra images, radially averages the correlation map, and identifies the nearest peak radius (characteristic spacing between Actin and Kibra signals).

**Requires:** Image Processing Toolbox, Signal Processing Toolbox

**Key parameters to set before running:**

| Parameter | Description |
|---|---|
| `Actin_img_name` | Actin surface projection image (TIFF) |
| `Kibra_img_name` | Kibra surface projection image (TIFF) |
| `Cellpose_img_name` | Cellpose segmentation mask (PNG) |
| `Output_file_name` | Output Excel file name |
| `Pixel_size` | Pixel size in µm/pixel |
| `Min_cell_area` | Minimum cell area in pixels (smaller cells excluded) |
| `Erode_iteration` | Number of erosion iterations to exclude junctions |

**Outputs:** Excel file with two sheets:
- `radial_mean_xcorr` — radially averaged cross-correlation per cell vs distance (µm)
- `local_peak_distance_after_mean` — nearest peak radius per cell (µm)

**Workflow:** Run `EVL_surface_projection.m` first to obtain the projected Actin and Kibra images, then segment cells with Cellpose, then run this script.
