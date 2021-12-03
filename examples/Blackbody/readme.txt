Here you will find the files of the Blackbody test done in Romero et al. (in prep), section 3.

-BlackBody_positions.txt contains the positions where the mean intensity output will be saved. This goes to the 'MeanIntensity_Positions' folder in 'input_data'
-BlackBody_source.stab is the emission spectra of the blackbody at the center of the simulation. This goes to 'star_sources' folder in 'input_data'
-BlackBody_gas.dat contains the initial conditions of the gas. The files Mixclask will read the input spectra are in the first column, from the folder than contains 'main.py'.
-Blackbody_star.dat contains the geometry and location of the stellar sources. Mixclask will look if the file in the first column exists in 'star_sources'.

To run this test, enter in 'Main.py' and write in these lines:
> # Where are the gas and star parameters starting from this folder?
> gas_params  = 'input_data/params/BlackBody_gas.dat'
> star_params = 'input_data/params/BlackBody_star.dat'
> meanIntensity_positions = 'input_data/MeanIntensity_Positions/BlackBody_positions.txt'
Assuming that you have placed these files as told above.

Now you can select in the "WavelengthOptions" dictionary your spectral resolution. For this test, I recommend to change 'convWavelength' from 150e3 to [5.0,150e3].
This option tells Mixclask what are the wavelengths where we test if the result obtained in a iteration converges.
For this case, We choose 5 nm and 150 microns because we are interested in having gas absoption and dust emission, respectively, correct. 

Finally, go to these lines:
> # Some technical parameters
> cloudy_path = '/path/to/your/cloudy/exe'
> show_cloudy_params = False
> last_iteration = 0
> n_iterations = 3
> n_threads = 4
In 'cloudy_path' you must specify where 'cloudy.exe' is found in your computer. 
If you installed cloudy normally, it should be in 'cXX.XX/source/cloudy.exe'. Where 'XX.XX' is the version of cloudy (e.g.: if version is 17.01, write c17.01).
PLEASE NOTE that this isn't the full path, you must put before the location of 'cXX.XX/source/cloudy.exe' (i.e.: the folder where you have installed cloudy), and the string must start with '/'
Also, change 'n_iterations' to a number higher than 10 (6 is okay for galaxies, but not for this test, I suggest to write n_iterations = 20 to follow Romero et al.), 
and change 'n_threads' to a lower value if you don't want to use a lot of CPUs.

And that's all, run in terminal:
> python3 Main.py
(just 'python' for newer versions, or run with spyder) and wait to the results.
This test can take some hours to complete, between 10 to 20h.
