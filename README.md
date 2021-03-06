# JawCloseR

JawCloseR is an assembly of two R scripts that, with a few other text files, models lower-jaw adduction in fishes using a dynamic equilibrium model based on that of van Wassenbergh et al. (2005), Kenaley (2012), and finally Kenaley et al. (2019). The "JawModel" script defines the function "run.jaw.prey.2." and runs the simulation interations with the "jaw.torque" file. These must be accompanied by a CSV input file of morphometric data and the parameter value file ("ParValPrey.csv"). For an example of the morphological data "C_sloani_mdata_for_MS.csv" and Kenaley et al. (2019) for explainations. The "var.df" file defines the variables that will be calculated during simulaitons.  The "ParValPrey" file defines imporatant simulation and muscle values, including the time step, muscle density, activation rise time, etc. The function "run.jaw.prey.2" input includes: 

* **spec.n**: Numeric, the specimen number identified in the "Specimen" column of the input file.
* **dat**: Character, the name of the CSV file.
* **config**: Numeric, the configuration of the adductors. "1" specifies no $A_{\omega}$ geometry, but its mass is added to PCSA calculations for the $A_2$. "2" specifies completely ignoring the $A_{\omega}$ and "3" specifies that the $A_{\omega}$ geometry and mass will be accounted for. These configuations should prove useful for many typical ray-finned muscle geometries.
* **Loose**: Logical, is it a loosejaw dragfish? Defaults to FALSE 
* **OutPutVerb**: Logical, should the simulations output all of the model calculations.
* **Prey**: Logical, should a prey item be added to the simulation. See below.
* **press**: Numeric, maximum intraoral pressure.
* **prey.per**: Numeric; how large (its length) is the prey item as a percetage of the specimen's body length.
* **strike.ang**: Numeric, when in adduction (the angle) the prey hits the jaw.
* **prey.pos**: Character, the position of the prey. Must either be "flat" (lateral side flat against jaw) or "flipped" (ventral or dorsal side against jaw.
* **MaxIts**: Numeric, how many interations of 0.1 MS should the simulation run.
* **progress**: Logical, should the progress be printed as a graphic stick figure of the muscle and jaw geometry. Printed automatically to a "\graphics" directory in working directory.
* **print.every**: Numeric, specifies how often to print progress.
* **out**: Character, output graphics file type. Must be either "png" or "pdf".

The model can accomadate the inclusion of a prey item against the jaw. This is very much in development and it would be best to ignore this option for the time being.
