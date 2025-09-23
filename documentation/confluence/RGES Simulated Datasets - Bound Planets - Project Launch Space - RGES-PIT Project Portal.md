#### Generated to Evaluate Science Filter Selection and Verify Mass Measurement Requirements

This page presents simulated microlensing datasets created to support ongoing studies of science filter selection and mass measurement requirements for the RGES. These simulations were produced under the _GBTDS Core Community Definition Committee_ \- recommended _"**overguide**"_ observing scenario, which uses:

-   _**Cadence:**_ 12.142 minutes
    
-   _**Exposure time:**_ 66.88 seconds
    
-   _**Fields:**_ 5 continuous fields + 1 Galactic Center field
    
-   **Color Cadence**: 6-hour
    

Each dataset assumes a 6-hour color cadence and realistic observational strategies using various combinations of Roman science filters. The simulations span a broad range of planetary masses and incorporate **Fisher matrix-derived uncertainties** for key microlensing parameters. To capture a representative range of planetary separations, **semimajor axes were drawn from a uniform distribution in logarithmic space between 0.3 and 30 AU**.

To ensure computational efficiency and focus on events with a non-negligible chance of planetary detection, we adopt the **caustic region of influence (CROIN) parameterization** framework described in _Penny (2014)_.

---

### **Planetary Mass Range**

To evaluate detection sensitivity and characterization performance across a broad range of planetary systems, we simulated events with the following planet masses:

| **Mass \[Earth Mass\]** | 0.1 M⊕ | 1M⊕ | 10M⊕ | 100M⊕ | 1000M⊕ | 10000M⊕ |
|---------------------|---------|-----|------|-------|--------|---------|

---

### Filter and Observatory Codes

The observatory codes used in these simulations correspond to Roman science filters:

| Observatory Code | Filter Name | Legacy Name |
|------------------|-------------|-------------|
|        0         |    F146     |    W146     |
|        1         |    F062     |    R062     |
|        2         |    F087     |    Z087     |
|        3         |    F184     |    F184     |
|        4         |    F213     |    K213     |
|        5         |    F106     |    Y106     |

---

## Observing Groups and Filter Combinations

All simulations were run for a 6-hour cadence. Each observing group corresponds to a unique filter trio:

| ObsGroup Number |   ObsGroup Name    | Filters \[F146 , blue , red\] |
|-----------------|--------------------|-------------------------------|
|        0        | overguide\_6hcc\_wrf |       F146, F062, F184        |
|        1        | overguide\_6hcc\_wrk |       F146, F062, F213        |
|        2        | overguide\_6hcc\_wzf |       F146, F087, F184        |
|        3        | overguide\_6hcc\_wzk |       F146, F087, F213        |
|        4        | overguide\_6hcc\_wyf |       F146, F106, F184        |
|        5        | overguide\_6hcc\_wyk |       F146, F106, F213        |

---

## Fisher Matrix Parameters: Linear vs Log Space

Fisher matrix uncertainties are computed for the following parameters, with their respective parameterization spaces:

| Parameter |    Space    |          Description          |
|-----------|-------------|-------------------------------|
|    t0     |   Linear    | Time of maximum magnification |
|    tE     | Logarithmic | Einstein radius crossing time |
|    u0     |   Linear    |       Impact parameter        |
|   alpha   |   Linear    |    Source trajectory angle    |
|     s     | Logarithmic |     Projected separation      |
|     q     | Logarithmic |          Mass ratio           |
|    rho    | Logarithmic |     Finite source effect      |
|   piEN    |   Linear    |  Parallax (North component)   |
|   piEE    |   Linear    |   Parallax (East component)   |
|    Fb     |   Linear    |         Baseline flux         |
|    fs     |   Linear    |          Source flux          |

_Note_: Flux parameters F<sub data-renderer-mark="true">b</sub> and f<sub data-renderer-mark="true">s</sub> are computed **individually for each filter**.

---

## Light Curve Sample

Each mass bin in the dataset includes a set of representative light curve samples. These events are a subset of those used in the **mass measurement requirements study**.

For visualization, each folder also includes a plotting file — `gulls_visual_diagnostics.pdf` — generated using the publicly available post-processing tool `gulls_viz`:  
🔗 [![](https://github.com/fluidicon.png)GitHub - gulls-microlensing/gulls-postprocessing: A repository for examples of gulls postprocessing and utility functions](https://github.com/gulls-microlensing/gulls-postprocessing/tree/main)

**Figure**: Example of bound planet light curve and caustic geometry from the `lc_sample/` subset. This event (6f\_overguide\_m10\_3\_5\_203.det.lc) is part of the 10 M⊕ sample. The **left panels** show the full and zoomed-in light curves featuring the planetary signal, comparing the best-fit planetary model (black) with a single-lens fit (red dashed). The **bottom-right panel** displays the Einstein ring (green dashed), caustic locations (purple dots), and the source trajectory (blue line), with the red arrow marking the source position at peak magnification (t0t\_0t0). The **top-right panels** show the detailed caustic structures.

---

## 🗂️ Folder Structure Notes

💡 Each top-level folder (6f\_overguide\_mXX/) corresponds to a single planet mass bin.

📂 The `analysis/` folder contains simulated data catalogs.

📦 The `det_lcout/det_lcout.tar.gz` file archives ~10% of all detected light curves.

📄 The `lc_sample/` folder contains a curated subset of sample light curves used for mass measurement requirement studies.

`filter_selection/ ├── 6f_overguide_m-10/ # Dataset for 0.1 M⊕ Bound Planet Events │ ├── analysis/ │ │ ├── 6f_overguide_m-10.out.hdf5 # Full event catalog │ │ ├── 6f_overguide_m-10.det.hdf5 # Detected events catalog (Δχ²>160) │ │ └── 6f_overguide_m-10.det.rates # Detection rates │ ├── det_lcout/ │ │ └── det_lcout.tar.gz # 10% of detected light curves │ ├── lc_sample/ │ │ ├── *.det.lc # Sample light curve files │ │ └── gulls_visual_diagnostics.pdf # Multi-page visualization using gulls_viz ├── 6f_overguide_m00/ # Dataset for 1 M⊕ Bound Planet Events ├── 6f_overguide_m10/ # Dataset for 10 M⊕ Bound Planet Events ├── 6f_overguide_m20/ # Dataset for 100 M⊕ Bound Planet Events ├── 6f_overguide_m30/ # Dataset for 1000 M⊕ Bound Planet Events ├── 6f_overguide_m40/ # Dataset for 10000 M⊕ Bound Planet Events`

---

## Data Access

You can access this dataset via [BOX](https://lsu.box.com/s/c9gvpr1h0fv8gls37sjonsfrcyunx7t5 "https://lsu.box.com/s/c9gvpr1h0fv8gls37sjonsfrcyunx7t5").