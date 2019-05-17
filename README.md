SDGPSR is a project to use a Software-Defined Radio and PC as a GPS Receiver. The codebase leverages *Understanding GPS: Principles and Applications, Kaplan and Hegarty* as well as IS-GPS-200J.

### Implementation
The receiver takes in data at 1 ms intervals (the length of the C/A sequence), and conducts the following steps in order:
1. Search the sky for all 32 GPS Satellites. A 2-dimensional search grid is created over code-phase and frequency offset, with a 10 KHz search window and 500 Hz increments. For each frequency, a circular correlation function is computed in the frequency domain, and the results are non-coherently integrated in 1 ms intervals across 128 ms.
2. A satellite is considered found if the value of one of the search bins is at least 10x the standard deviation of all of the search bins. 
3. For each satellite that is found, a tracking channel is created. Each tracking channel is initialized with multiple signal trackers, all with different frequency offsets. 
  * Each signal tracker attempts to lock onto the signal
  * Once one signal tracker obtains the full track state, the others are removed
4. The signal tracker proceeds through 4 states in order to get a signal lock
  * It starts with 1 ms integration intervals and attempts to get a frequency lock on the carrier
  * Once a frequency lock is completed, it attempts to get a phase lock on the carrier
  * Once carrier phase lock is completed, it starts looking for nav bit edges. These are where the carrier changes phase by 180 degrees.
  * Once the nav bit edges are found, it starts demodulating the nav bits
5. The LNAV data is demodulated. 
  * In order to find where we are in a subframe, we search for the preamble on word 1. If the preamble is found (or its inverse, indicating we have locked phase 180 off), an attempt to validate the word with its checksum is made. 
  * If this is successfull, a new subframe demodulation is started. 
  * If the nav bits are found to be flipped, a flag is returned to flip the phase on the signal tracker.
  * Time of week is determined from any subframe
  * Subframes 1-3 are demodulated, and the values populated in structures for clock data and ephemeris data. Subframes 4 & 5, which contain the Almanac and other less essential information are ignored
6. The clock and ephemeris data is used to determine the position of the satellite in the ECEF frame
7. The satellite positions and time of arrival of each signal is used to compute the user time & position using least squares

### Major items that still need to be implemented:

* Satellite Errors (clock, relativistic, ionospheric delays, tropospheric delays, etc)
* Sagnac effect
* Add/Remove satellites as they come into/disappear from view
* Demodulate Almanac to determine when new satellites are visible
* Satellite velocity
* Mapping Doppler Measurements to Satellite Velocity for user velocity
* Implementation of the Wide Area Augmentation System (WAAS)
* Code cleanup from original development
* Makefile for shared library
* Writeup of the process and results
