# OH-Calibration-for-HONO-Cell
Author: Bode Hoover (bodehoov@iu.edu)
Adapted from VBA code written by Emily Reidy
Last updated: July 5, 2022

Performs OH calibration to determine sensitivity parameters. 
Inputs laser power, relative humidity, temperature, online &amp; offline signal, and other variables. Calculates OH and HONO LOD.




Default extractor constants must be edited prior to running (laser powers and UVouts). 
Must be txt file where spaces are replaced with commas and comma added after MeasCycle column.
Must change directory for unique user for output.
Calibrator characterization constants may be updated periodically. 

User is prompted for 2 inputs:
1. Calibrator #
2. Water source (box monitor or probe)

Script determines cycles numbers, averages each cycle, separates online & offline, calculates both linear and non linear best fits, calculates OH sensitivity, calculates OH and HONO LODs.
