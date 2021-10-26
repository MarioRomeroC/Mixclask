These are the folders used by mixclask.

-gas_props and gas_sources saves the opacity and emission spectra, respectively, of the gas. Their files are generated and used by this code, without user intervention.
-star_sources stores the emission spectra of any stellar sources, and should be included by the user. All files there should end in '.stab' (e.g.: your_source.stab)
-params includes the input data of the run, and MeanIntensity_positions the positions of your simulation where you want the mean intensity. Info for these files are not given in this readme.

About files in star_sources, they have this header:
 # Column 1: wavelength (nm) 
 # Column 2: specific luminosity (erg/s) 
 # Normalization wavelength X
 # Value Y
Where X is a wavelength in nm, and Y is the value of the specific luminosity (in erg/s) in that wavelength. With specific luminosity I mean nu*L_nu = lambda*L_lambda.

Files inside MeanIntensity_positions should have this structure:
(1) The header is:
 # Positions to export the mean intensity field
 # column 1: position x (pc)
 # column 2: position y (pc)
 # column 3: position z (pc)
This is the same file you import in skirt to get the probe of radiation field at positions (https://skirt.ugent.be/skirt9/class_radiation_field_at_positions_probe.html) and the same rules apply here (for example, you can change 'pc' by 'kpc' to give the column in those units).
(2) However, the order of the rows are very important here. You must give first the positions of the zones you are modelling. For example, if you model a galaxy as rings centered in R=1,3,5 kpc (and z=0), first three rows are:
 1000.0 0.0 0.0
 3000.0 0.0 0.0
 5000.0 0.0 0.0
If you use the same header as above. These files will be saved with the same name as the input sed for the previous iteration (or your initial condition), stored in the input data given for the gas, inside the params folder.
(3) Next rows are any custom positions you want, if desired. They will be saved as 'extra_[r]pc.sed', where '[r]' is sqrt(x^2 + y^2 + z^2)


