# Native PyMOL: three DENSS representations

These are actual PyMOL renders of WT_refined.mrc and AA_refined.mrc from the DENSS reconstruction/refinement workflow. Both samples use the same camera scale and two views separated by 90 degrees.

- Gray space-filling: density beads on a common 4 Å grid, sphere radius 2 Å, inside the contour enclosing 90% of positive density. These are not atoms or an atomic model.
- Blue: translucent isosurface enclosing 90% of positive density.
- Rainbow: continuous OpenGL volume rendering, using the same ramp expressed as fractions of each map's maximum. Blue is lower and red higher relative density. Colors are not absolute-density comparisons between samples.

One PyMOL map_double interpolation is used for display. No structural model is overlaid. The maps themselves are unchanged. For this IDP, the maps are exploratory effective scattering envelopes, not a unique conformation; lobes are not established domains.

Open DENSS_three_styles.pze in PyMOL and select the named scenes. To reproduce, run setup_views.py followed by export_views.py in the PyMOL GUI; rainbow requires an OpenGL context. assemble.py only lays out the resulting images. See render_settings.json for exact thresholds and ramps. Native software: DENSS 1.8.7, PyMOL 3.1.0.
