# Mixclask
Mixing cloudy and skirt

This code combine cloudy (https://gitlab.nublado.org/cloudy/cloudy/-/wikis/home) and skirt (https://skirt.ugent.be/root/_landing.html) to predict spectra and gas properties in astrophysical contexts (galaxies, HII regions...). This code assumes that you have installed both codes following their instructions.
You must also have python3 installed with flatten-dict module (if not, write 'pip3 install flatten-dict' in terminal)

The main output of this code is the mean intensity
> 4π*λ J = λ int( I dΩ)
of a region filled with stars, gas and dust; at different positions (R,z), assuming axial symmetry.
Mixclask divides the simulation into different stellar and ISM (i.e.: gas+dust) regions, with different properties.


In order to run Mixclask, you need as inputs the stellar and ISM data for each region, and an file indicating the positions (x,y,z) where you want to generate an output.
These files have an structure of header+data, such as
># column 1 : param1 (units1)
># column 2 : param2 (units2)
># column 3 : param3 (units2)
>...
>data11 data12 data13 ...
>data21 data22 data23 ...
>...
dataXX are usually some number, but it can be a python-list or a filename as well. Different data in a row is ALWAYS separated by a space.

For particular input files, you have.
-Stellar info:
You need to provide two columns per row: an spectra for that region, and its geometry.
The spectra is another file containing the neutral luminosity λ*L and wavelength λ. Following this structure:
># Column 1: Wavelength (nm) 
># Column 2: Luminosity (erg/s) 
># Normalization wavelength 550.0 
># Value 1.17e+43 
>9.1 2.26e+37 
>...
As you see, you need two additional header lines. First one is what it says, in nm. Second line is the value of λ*L (erg/s) at normalization wavelength (in the example λ*L=1.17e+43 erg/s at 550 nm).

With the stellar spectra files, the stellar input data is, for example:
># Column 1 : SedFile 
># Column 2 : Geometry ['type',params](pc)
>file1.stab [’point’,0.0,0.0,0.0]
>file2.stab [’ring’,1.0,0.1,0.2]
>file3.stab [’shell’,2.0,3.0] 
Mixclask consider each row a different region with different data. SedFile is the spectra filename, LOCATED INSIDE 'input_data/star_sources'.
This direction can be changed in skirt/ski_params.py file, in the variable 'self._star_sources_folder' .

Geometry is a python-list containing their geometry and its parameters (separated by commas).
For stars, you have three options.
1-Point geometry, the remaining parameters are the (x,y,z) position of the region. Right now, it only accepts points at the origin (sorry)
2-Ring geometry, assumes a ring centered in the origin with radius, width and height as remaining parameters in the list (see https://skirt.ugent.be/skirt9/class_ring_geometry.html for details)
3-Shell geometry, assumes an spherical shell with constant values centered in the origin, with inner radius and outer radius as remaining parameters (https://skirt.ugent.be/skirt9/class_shell_geometry.html for details)
In the example I provided, the simulation domain have three stellar regions: First one is a point source located in the origin. Second one is a ring of radius 1.0 pc, width 0.1 pc and height 0.2 pc. And third region is a shell between 2.0 and 3.0 pc.

-ISM info:
Similar structure as stellar info, but you have several and optional parameters to add.
An example for this input file may be
># Column 1 : SedFile 
># Column 2 : Mass (Msun) 
># Column 3 : HydrogenDensity (cm-3) 
># Column 4 : Geometry ['type',params](pc) 
># Column 5 : Helium fraction (1) 
># Column 6 : Metallicity (1) 
># Column 7 : Dust-to-gas ratio (1) 
># Column 8 : PAH-to-dust ratio (1) 
>ISMregion1.sed 1e5 1.0 ['ring',1.0,0.5,0.2] 0.24 0.02 0.02 0.05
>ISMregion2.sed 1e7 1.0 ['shell',2.0,3.0] 0.24 0.02 0.02 0.05
The first four columns are mandatory to run Mixclask. 
SedFile here is the name of the spectra that Mixclask generates for each ISM region (you don't need to provide an extra file, is the output filename of that region).
Mass is the total ISM mass (gas+dust) of each region, and HydrogenDensity (without spaces) the number of Hydrogen particles per cm^-3.
Geometry follows the same structure as the stellar geometry for ring and shell geometries (point geometry is disabled for ISM regions).

For chemical elements, you can give more data such as He fraction and Z. These are mass fractions relative to ONLY GAS.
Furthermore, you can give specific elements mass fractions up to Zn. They are labelled as 'lithium fraction', 'beryllium fraction', 'boro fraction', 'carbon fraction', etc.
If you give specific elements, you must not give metallicity, as it takes precedence.
By default, Mixclask uses Grevesse et al. 2010 abundances (Solar composition, 'abundances gass' in cloudy) as base.

For dust, you can give the Dust-To-Gas ratio, M_dust/ISM mass and the PAH-to-dust ratio (M_PAH/M_dust).
PAH-to-dust ratio is named in the literature as q_pah.
By default, Mixclask adds ism grains and PAH, using the cloudy commands 'grains ism' and 'grains pah' (see cloudy documentation). To remove dust, you have to give DTG = 0.0 in the input.

-Output positions: 
File structure is the same as indicated in 'AtPositionsForm': https://skirt.ugent.be/skirt9/class_at_positions_form.html
># column 1: position x (pc)
># column 2: position y (pc)
># column 3: position z (pc)
>1.0 0.0 0.0
>2.5 0.0 0.0
>0.0 0.0 1.0
>3.0 2.0 5.0
>...
Units of each columns are more flexible here (you can write 'kpc' instead of 'pc' and give data in these units)
However, first rows ALWAYS have to match with the middle of each ISM region. For rings, this is the ring radius. For shells, this is the sum of the outer and inner radius, divided by two (its mean).
In the ISM example provided, first row must be 1.0 0.0 0.0 and second 2.5 0.0 0.0. The output files of these region would be the SedFile you provide in the ISM input data ('ISMregion1.sed' and 'ISMregion2.sed', from the same example).
Once all ISM regions have its corresponding row, any remaining rows can have the positions you want, and they will be named as 'output_R[number]_z[number]kpc.sed .
Taking the ISM input example, third row will produce an output called 'output_R0.0036_z0.005kpc.sed'

---

Once you have all input files, next you have to enter in 'Main.py'.
There, you find a dictionary called 'Options', subdivided in several dictionaries representing a category.
First, you need to configure these lines:
>'FileParameters':{
>    'gas_params':'input_data/params/Your_gas_file.dat',
>    'star_params':'input_data/params/Your_star_file.dat',  
>    'positions':'input_data/Output_Positions/Your_positions.txt'
>},
These contains the location of the input data files explained above.

Next, you set the wavelength resolution:
>'Wavelength':{ #Wavelength related options.
>    #Resolution and wavelength range are specified here
>    'maxWavelength':3.0e5, #nm
>    'minWavelength':10.0, #nm
>    'resolution':200, #Equispaced in log
>    #Other wavelength options
>    'normalization': 550.0, #nm -> This is used to normalize Skirt sources luminosity for ISM sources
>},
The wavelength range is defined with 'maxWavelength' and 'minWavelength', and the number of points between them with 'resolution'. Points are equispaced logarithmically.
As other options, 'normalization' is the normalization wavelength of the ISM emission spectra files that Mixclask generates. Although any number between 'maxWavelength' and 'minWavelength' works, I recommend to use the same wavelength used for stellar spectra files for consistency.

Then, you have options related with convergence and speed
>'AccuracyAndSpeed':{
>    #Convergence options
>    'convWavelength':[(10.0,90.0),(100e3,300e3)], #nm #Wavelength (value/range) to look convergence
>    'tolerance': [0.67,0.10], #(WIP) How much do you want to differ from previous iteration results. Must have the same length as 'convWavelength'
>    #Speed options
>    'n_threads': 2, #Number of logical cores you want to run for a SINGLE simulation
>    'photon_packets':1e8
>},
First, 'convWavelength' is the set of wavelengths, or wavelength ranges, in which Mixclask will check for convergence.
You can provide a list such as, for example:
>'convWavelength':[(10,90),300,(500,600),150e3]
With this list, Mixclask will check convergence in:
1-the range between 10 to 90 nm
2-300 nm
3-Range between 500 to 600 nm
4-150e3 nm
By convengence, Mixclask understands that the mean intensity should not vary between current and previous iterations for EACH ISM REGION AND EACH ELEMENT OF THE LIST (one or two elements in the list may be enough).
That is, the value does not change when you provide a single wavelength (300 and 150e3 nm, in this case), or the integral in the range provided does not change (for (10,90) and (500,600) elements).
The maximum amount allowed to change is given in 'tolerance'. 
For example, if 'convWavelength' = [(10.0,90.0),(100e3,300e3)] and tolerance = [0.67,0.10], Mixclask will check tolerance = 0.67 for the (10.0,90.0) nm range, and 0.1 for the (100e3,300e3) nm . 
Convergence is evaluated with this expression: 
> |(This iteration) - (Prev. iteration)| <= tolerance*[(This iteration) + (Prev. iteration)]/2 
for all ISM regions and wavelengths given in the ISM file and 'convWavelength', respectively
On the other hand, 'n_threads' and 'photon_packets' are parameters needed for skirt. The former is the number of task that your cpu will use for executing this code (usually, n_threads = 2*n_cpu).
The latter is the number of photons used in the Montecarlo simulation of skirt, the bigger the number, less noiser the results (but also slower).

Finally, you also have 'Technical' options. The only one you have to touch is
> 'cloudy_path':'/path/to/your/cloudy/exe' 
In 'cloudy_path', you must specify where 'cloudy.exe' is found in your computer. 
If you installed cloudy normally, it should be in 'cXX.XX/source/cloudy.exe'. Where 'XX.XX' is the version of cloudy (e.g.: if version is 17.03, write c17.03).
PLEASE NOTE that this isn't the full path, you must put before the location of 'cXX.XX/source/cloudy.exe' (i.e.: the folder where you have installed cloudy), and the string must start with '/'.

---

And that's all! To run Mixclask, just write in terminal
>python3 Main.py
Also, you can run it in spyder instead if you wish.
As an alternative, you can run it in background as
> python3 -u Main.py > file.log &
