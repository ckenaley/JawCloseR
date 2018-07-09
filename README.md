# JawCloseR

JawCloseR is an assembly of two R scripts that, with a few other text files, models lower-jaw adduction in fishes using a dynamic equilibrium model based on that of van Wassenbergh et al. (2005) and Kenaley (2012). The "jaw.model" script defines the function "run.jaw.prey.2." It must be accompanied by a CSV input file of morphometric data. See example under "C_sloani_mdata_for_MS.csv" and Kenaley et al. (submitted) for explainations. Must also include  The function's input includes: 

* **spec.n**: Numeric, the specimen number identified in the "Specimen" column of the input file.
* **dat**: Character, the name of the CSV file.
* **config**: Numeric, the configuration of the adductors. "1" specifies no $A_{\omega}$ geometry, but its mass is added to PCSA calculations for the $A_2$. "2" specifies completely ignoring the $A_{\omega}$ and "3" specifies that the $A_{\omega}$ geometry and mass will be accounted for. These configuations should prove useful for many typical ray-finned muscle geometries.
* **Loose**: Logical, is it a loosejaw dragfish? Defaults to FALSE 
* **OutPutVerb**: Logical, should the simulations output all of the model calculations.
* **Prey**: Logical, should a prey item be added to the simulation.
* **press**: Numeric, maximum intraoral pressure.
* **prey.per**: Numeric; how larger is the prey item as a percetage of the specimens body length.
* **strike.ang**: Numeric, when in adduction (the angle) the prey hits the jaw.
* **prey.pos**: Character, the position of the prey. Must either be "flat" (lateral side flat against jaw) or "flipped" (ventral or dorsal side against jaw.
* **MaxIts**: Numeric, how many interations of 0.1 MS should the simulation run.
* **progress**: Logical, should the progress be printed as a graphic stick figure of the muscle and jaw geometry. Printed automatically to a "\graphics" directory in working directory.
* **print.every**: Numeric, specifies how often to print progress.
* **out**: Character, output graphics file type. Must be either "png" or "pdf".
