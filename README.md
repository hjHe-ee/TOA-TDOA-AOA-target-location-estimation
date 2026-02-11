# TOA-TDOA-AOA-target-location-estimation

This repository provides **MATLAB-runnable examples** of three classic 2D target localization methods:

- **TOA** (Time of Arrival/Equivalent Distance)

- **TDOA** (Time Difference of Arrival/Hyperbolic Localization)

- **AOA** (Angle of Arrival/Direction Finding)

Each method constructs a simulation scene under **minimum number of base stations** conditions, generates noisy measurements, and performs:

1) Position estimation (True vs Estimated)

2) Visualization of geometric constraints (dashed lines: circle / hyperbola / azimuth ray)

3) Command-line output of estimation results and errors

> Note: The geometric trajectory of TDOA in a 2D scene is a **hyperbola (constant distance difference)**, not an ellipse.
