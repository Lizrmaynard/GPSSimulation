# GPSSimulation
Final project for numerical analysis MATH 5610.

This is a set of codes to aid in simulating a satellite array in which a pathway from one point to another can be mapped on the globe. In satellite.py, the code takes data from 24 satellites and determines which four satellites are optimal to use Newton's method to find the location of the individual, and continually what of the 24 satellites are appropriate to use for each step of the mapped pathway.

receiver.py then takes the satellite data given by satellite.py (via a provided vehicle by the instructor) and calculates the locations of each step and gives an output of coordinate locations on a globe for each step.
