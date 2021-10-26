Here you will find the files of the Blackbody test done in Romero et al. (in prep), section 3.

-BlackBody_positions.txt contains the positions where the mean intensity output will be saved. This goes to the 'MeanIntensity_Positions' folder in 'input_data'
-BlackBody_source.stab is the emission spectra of the blackbody at the center of the simulation. This goes to 'star_sources' folder in 'input_data'
-BlackBody_gas.dat contains the initial conditions of the gas. The files Mixclask will read the input spectra are in the first column, from the folder than contains 'main.py'.
-Blackbody_star.dat contains the geometry and location of the stellar sources. Mixclask will look if the file in the first column exists in 'star_sources'.
-Iteration0 contains the input files that are described in 'BlackBody_gas.dat'.

To run this test, COPY all files from Iteration0 into the root folder (where is 'Main.py').
Enter in 'Main.py' and write in these lines:
> # Where are the gas and star parameters starting from this folder?
> gas_params  = 'input_data/params/BlackBody_gas.dat'
> star_params = 'input_data/params/BlackBody_star.dat'
> meanIntensity_positions = 'input_data/MeanIntensity_Positions/BlackBody_positions.txt'
These assumes that you have placed these files as told you above.

Now you can select in the "WavelengthOptions" dictionary your spectral resolution. For this test, you don't need to change these lines.
Finally, go to these lines:
> # Some technical parameters
> cloudy_path = '/path/to/your/cloudy/exe'
> show_cloudy_params = False
> last_iteration = 0
> n_iterations = 3
> n_threads = 4
In 'cloudy_path' you must specify where 'cloudy.exe' is found in your computer. 
If you installed cloudy normally, it should be in 'cXX.XX/source/cloudy.exe'. Where 'XX.XX' is the version of cloudy (e.g.: if version is 17.01, write c17.01).
PLEASE NOTE that this isn't the full path, you must put before the location of 'cXX.XX/source/cloudy.exe' (i.e.: the folder where you have installed cloudy)
Also, change 'n_iterations' to a number higher than 10 (3 is okay for galaxies, but not for this test, I suggest to write n_iterations = 15 to follow Romero et al.), 
and change 'n_threads' to a lower value if you don't want to use a lot of CPUs.

And that's all, run in terminal:
> python3 Main.py
(just 'python' for newer versions, or run with spyder) and wait to the results.
This test in my laptop took about 6-10 hours.
