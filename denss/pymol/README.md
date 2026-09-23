# DENSS maps rendered with PyMOL

Rendered by PyMOL Open-Source 3.1.0 using its native isosurface and ray-tracing commands. The workflow follows DENSS's visualization guide: https://tdgrant.com/tips/ . This replaces the previous custom-renderer presentation.

The PNG/PDF figures show three translucent isodensity surfaces: blue outer contour (90% of positive density), cyan (60%), and red core (30%). These fractions define display thresholds independently for WT and AA, not statistical confidence or absolute density differences. Identical camera settings are used for both samples. The density maps are the existing refined maps; no reconstruction or fit was rerun. PyMOL map_double interpolates a display copy of each map; the source MRC files are unchanged.

`DENSS_volume_ramps.pze (select scene `raytraced_shells`)` opens the native ray-traced shell scene. `DENSS_volume_ramps.pze` opens continuous volume coloring, using the volume/volume_ramp_new approach described by DENSS. The continuous volume display requires PyMOL's OpenGL viewer; it is not the representation in the ray-traced PNGs. The saved scenes allow rotation and threshold adjustment. The colored cores are density features of an effective reconstruction, not identified protein domains.

Run `render_pymol.py` using a Python environment with pymol-open-source and numpy to regenerate the native images and sessions. Run `assemble.py` (Pillow and matplotlib) for figure layout. No custom renderer is used for the maps; the layout script only places the images and labels.
