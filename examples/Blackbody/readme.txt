Here you will find the files of the Blackbody test done in Romero et al. (2023), section 2.3.

-BlackBody_positions.txt contains the positions where the mean intensity output will be saved.
-BlackBody_gas.dat contains the initial conditions of the gas. The files Mixclask will read the input spectra are in the first column, from the folder than contains 'main.py'.
-BlackBody_star.dat contains the geometry and location of the stellar sources. Mixclask will look if the file in the first column exists in 'star_sources'.
-Inside the 'star' folder, you will find the stellar spectra of all regions in 'BlackBody_stars.dat'.
You may move the contents to 'input_data' folder, but that's not longer necessary.

To run this test without minimal input for the user, enter in 'Main.py' and write in these lines:
>'FileParameters':{
>        'stars':{
>            # Where is the star parameter file?
>            'file': 'examples/Blackbody/BlackBody_star.dat',
>            # And the folder containing the SED of each stellar region?
>           'folder': 'examples/Blackbody/star'
>        },
>        # Where is the ISM parameter file?
>        'ISM': 'examples/Blackbody/BlackBody_gas.dat',
>        # At what positions (x,y,z) do you want some outputs?
>        'positions': 'examples/Blackbody/BlackBody_positions.txt'
>    },
Assuming that you have not moved these files.

Next, go to this dictionary and write:
>'PhotonProbability':{
>    'per_region':'Custom', # Available options:'logWavelength','Custom'
>        # This option modifies, inside each region, the probability distribution of which a photon of certain wavelength is launched
>        # Available options are:
>        # 'logWavelength': p(λ) ~ 1/λ -> Follows the Logarithmic distribution option given in Skirt
>        # 'Custom': p(λ) ~ f(λ) -> Given by the user below (used for ALL regions)
>    'customDistributionFile': 'input_data/probability_distributions/EUV_50percent.stab',
>    'wavelengthBias': 0.5 #Between 0 and 1. It controls how many photons launched per region follow above distribution.
>        # 0 makes above option without effect.
>        # 1 makes regions to strictly follow the distribution
>        #   (you risk having some wavelength ranges without photons, so the output will be zero there,
>        #       because the probability was too low to even launch one photon)
>}
I am changing 'per_region' and 'customDistributionFile'. Do not worry about the file, because it is inside 'input_data/probability_distributions'.

Finally, go to this line:
> 'cloudy_path':'/path/to/your/cloudy/exe',
In 'cloudy_path' you must specify where 'cloudy.exe' is found in your computer. 
If you installed cloudy normally, it should be in 'cXX.XX/source/cloudy.exe'. Where 'XX.XX' is the version of cloudy (e.g.: if version is 17.03, write c17.03).
PLEASE NOTE that this isn't the full path, you must put before the location of 'cXX.XX/source/cloudy.exe' (i.e.: the folder where you have installed cloudy), and the string must start with '/'.

Feel free to change other parameters such as:
> 'n_iterations' (maximum number of iterations if convergence is not reached earlier)
> 'n_threads' (number of logical cores you want to use in Skirt)
> 'n_cpus' (Number of simulations to be run at once in Cloudy)
> 'photon_packets' (number of photon launched in each skirt run)
These affect mainly computation speed. 

And that's all, run in terminal:
> python3 Main.py
(just 'python' for newer versions, or run it with spyder) and wait to the results.
If you want the text that Mixclask says during runtime written in a file instead, use
> python3 -u Main.py > file.log &

This test can take some hours to complete, about 7-10h.

---

Once it is finished, I suggest to create a folder named 'BlackbodyResults' (for instance, you may use other name) and move all folders called 'iteration' and any '.sed' files that are in the root folder.
If you want the average of all iterations, go to 'post_processing/averages.py'. Inside the file, change this line:
> inputs_path = '../BlackbodyResults/'
And then run
> python3 averages.py
Inside 'BlackbodyResults' folder, there should be another one called 'Statistics' with the average results.

