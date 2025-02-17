# Matilda
Live data processing package for USAXS/SAXS/WAXS instrument.

Not even started properly, yet.

**Planned features:**
For Bonse-Hart USAXS (https://usaxs.xray.aps.anl.gov) be able to reduce data collected using both Flyscanning and Stepscanning. For SAXS and WAXS using area detectors do the same while applying default masks. 
Reduce and normalize these data, apply thickness and absolute intensity correction (USAXS) or absolute intensity scaling parameter.

*If user provides proper name of "blank" (AKA: empty/blank/background) data BEFORE data collection (in instrumental epics interface), subtract the blank and generate data on absolute intensity scale - slit smeared and desmeared using Lake method.*



