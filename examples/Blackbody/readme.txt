Here you will find the files of the Blackbody test done in Romero et al. (in prep), section 3.

-BlackBody_positions.txt contains the positions where the mean intensity output will be saved. This goes to the 'MeanIntensity_Positions' folder in 'input_data'
-BlackBody_source.stab is the emission spectra of the blackbody at the center of the simulation. This goes to 'star_sources' folder in 'input_data'
-BlackBody_gas.dat contains the initial conditions of the gas. The files Mixclask will read the input spectra are in the first column, from the folder than contains 'main.py'.
-Blackbody_star.dat contains the geometry and location of the stellar sources. Mixclask will look if the file in the first column exists in 'star_sources'.

To run this test without minimal input for the user, enter in 'Main.py' and write in these lines:
> # Where are the gas and star parameters starting from this folder?
> gas_params  = 'input_data/params/BlackBody_gas.dat'
> star_params = 'input_data/params/BlackBody_star.dat'
> meanIntensity_positions = 'input_data/MeanIntensity_Positions/BlackBody_positions.txt'
Assuming that you have placed these files as told above.

Finally, go to these lines:
> # Some technical parameters
> cloudy_path = '/path/to/your/cloudy/exe'
> show_cloudy_params = False
> n_iterations = 15
> n_threads = 2
> tolerance = 0.1
In 'cloudy_path' you must specify where 'cloudy.exe' is found in your computer. 
If you installed cloudy normally, it should be in 'cXX.XX/source/cloudy.exe'. Where 'XX.XX' is the version of cloudy (e.g.: if version is 17.01, write c17.01).
PLEASE NOTE that this isn't the full path, you must put before the location of 'cXX.XX/source/cloudy.exe' (i.e.: the folder where you have installed cloudy), and the string must start with '/'.
Feel free to change the n_iterations (=maximum number of iteration if convergence is not reached earlier) and n_threads (~number of cpus you want to use) if you want. 

And that's all, run in terminal:
> python3 Main.py
(just 'python' for newer versions, or run it with spyder) and wait to the results.
If you want the text that Mixclask says during runtime written in a file instead, use
> python3 -u Main.py > file.log &

This test can take some hours to complete, between 10 to 20h.
