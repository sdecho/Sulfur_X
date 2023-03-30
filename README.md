# Sulfur_X
A model of sulfur degassing during magma ascent

Please cite Sulfur_X as follows:
Ding, S., Plank, T., Wallace, P., Rasmussen, D. J., in press. Sulfur_X: A model of sulfur degassing during magma ascent. Geochemistry, Geophysics, Geosystems. 
https://doi.org/10.31223/X56H0F. The related manuscritp is available on EarthArxiv: https://eartharxiv.org/repository/view/3559/

1.	Instructions for using Sulfur_X
Sulfur_X is written in Python and can be run in python 3.9 or any version higher. It requires packages of Pandas, Numpy and Scipy. The package includes one main function script, several other scripts of different classes called by the main function, and two example csv files for inputting melt inclusion data for comparison. All of these files should be kept in the same folder, and none of the names should be changed in any python script file. 

To run the program, open the main function script named “main_Fuego.py”, change the basic inputs and the advanced parameters accordingly and run. Melt composition can be modified in the python script named “melt_composition.py”. Parameters for low pressure (<20 MPa) degassing (INT and BAR) can be changed in the python script named “degassingrun.py”. The example input parameters that will reproduce the results for Fuego and Mauna Kea in the main text of this paper are listed in Table 1. 
1.1.	Basic input
temperature: Magma temperature in °C. Degassing is modeled isothermally. 

delta_FMQ: fO2 relative to FMQ buffer (Frost, 1991). If the redox evolution is enabled, this is the initial fO2 at the initial pressure. If the redox evolution is disabled, the degassing run will be run along this buffer. The initial Fe3+/FeT in the melt is calculated based on this input (Kress and Carmichael, 1991)

H2O_initial: Initial water concentration in the melt in wt.%.

CO2_initial: Initial carbon dioxide concentration in the melt in ppm.

S_initial: Initial total sulfur concentration in the melt in ppm.

choice: choice of crystallization or not. If equals 1, crystallization is enabled; if equals 0, crystallization is disabled. 

COH_model: choice of COH degassing model. If equals 1, VolatileCalc (Newman and Lowenstern, 2002); if equals 0, Iacono-Marziano model (Iacono-Marziano et al., 2012). We recommend using VolatileCalc for low water (<1wt.%) system.

fO2_tracker: choice of redox change along S degassing or not. If equals 1, fO2 changes by electron exchange between S and Fe in the melt and the co-existing vapor; if equals 0, fO2 is always buffered at the initial delta FMQ.

monte_carlo: Optional mote-carlo simulation for S degassing error estimate. If equals 1, monte carlo simulation will be performed after the first degassing run, and input the number of simulations as “m_run”; if equals 0, the program finishes after one degassing run. monte_name: Name (input as str and ends in .csv) of the output csv file of S concentration in the melt (in ppm) along decompression of each simulation, the mean, standard deviation and variance of S (in ppm) at each pressure from multiple simulations. This file will be rewritten each time the code is run (even if the monte-carlo simulation is turned off). In order to save the results, either change the name or save the output file under a different name before the next run. 
l: The total number of pressure steps between the initial pressure and 1bar. We recommend this number to be between 300-500. Too large or too small pressure steps can both cause unexpected excursions in the model results.
m: The total number of the degassing runs during decompression. 0 ≤ m ≤ l.

sulfide: Sulfide composition, S, Fe, Ni, Cu and O, in wt.%. This is used to calculate the sulfur contents at sulfide saturation (SCSS) using Smythe et al. (2017) model. The SCSS calculated here is only for the reference. Sulfur_X currently does not consider sulfide saturation along degassing.

1.2.	 Advanced parameters
Changes in any of the following parameters could have significant impacts on the model results. Please change the following parameters with caution and note the changes if the results are used.

1.2.1 key parameters for redox and crystallization calculation in main_Fuego.py

S_Fe_choice: choice of S speciation model. If equals 0, Nash model is employed (Nash et al., 2019); if equals 1, O’Neill and Mavrogenes (2022) model is employed; if equals 100, Muth model is employed (Muth and Wallace, 2021); if equals any other float number, revised Muth model is employed and the input float number is the last constant in the Muth model. The choice of this model is the key for modeling sulfur degassing and fO2 evolution. We recommend users choose the model that fits the S-Fe speciation measurements in the melt for the volcanic system of interest.

Sigma: The tolerance of the calculated log10fO2. When change in log10fO2 is smaller than sigma, the iterations searching for fO2 satisfying all the conditions stops (more details in the text). For the example runs, sigma equals to 0.05. If sigma is too small or too large, it may cause large excursions or even false fO2 calculations in the run.

If the crystallization is enabled, H2O (in wt.%)-melt fraction relation is specified using H2O (wt.%)-K2O (wt.%) correlation (K2O = a*H2O +b), assuming K2O is perfectly incompatible. The following two parameters, slop_h2o and constant_h2o represent a and b, respectively. 
slop_h2o: slop of H2O (wt.%) in K2O (wt.%) as a function of H2O.
constant_h2o: constant in K2O (wt.%) as a function of H2O. 
Change the slope and constant based on the H2O-K2O relation of the volcanic system of interest. If the system is assumed to be similar to Fuego, or crystallization is disabled, leave them unchanged. 

1.2.2. Melt composition in melt_composition.py
Melt composition is used to calculate the sulfur partition coefficients and COH degassing if the Iacono-Marziano model is used. Initial melt composition and the change as a function of melt fraction can be changed in the python script named “melt_composition.py”. If crystallization is enabled, the default input of the melt composition change as a function of melt fraction is based on Fuego magma (Lloyd et al., 2013; Rasmussen et al., 2020). If crystallization is disabled, the default input of the melt composition is an averaged composition for Hawaiian magma (Moussallam et al., 2016; Brounce et al., 2017). 
1.2.3. Low pressure degassing (<25 MPa) in “degassingrun.py”
Since there is very little experimental information at pressures lower than 25 MPa except for the H-free experiments from Nash et al. (2019) and O’Neill (2002), the predicted sulfur partition coefficients increase exponentially at this very low-pressure range and large uncertainties are involved. This might cause large excursions of the model results at low pressure, especially in the fO2 calculation. Therefore, Sulfur-X has the option to impose an arbitrary incremental increase in the combined molar partition coefficients of sulfur (INT) at each pressure step at low pressure (BAR in MPa, should be smaller than 20 MPa). INT and BAR can be changed in the python file named “degassingrun.py”. If there are large excursions in the fO2 calculation in the last a few pressure steps, we recommend the users to play with these two parameters to get the best fO2 results.
1.3 Input of melt inclusion data for comparison in main_Fuego.py
If you’d like to compare your model results to H2O-CO2-S-Fe3+/FeT-S6+/ST data measured in the melt inclusions, you can use either of provided csv files (named “Fuego.csv” and “Hawaii.csv”) as a templet and change the read-in csv name in main_Fuego.py. Please do not change any of the column names in the csv templet. 
mi_name: name (input as str) of the csv file containing melt inclusion data. 
1.4 Output file and plots in main_Fuego.py
All the results in the dataframe “df_results” are output in a csv file. The output file name should be changed for each run. Otherwise, the result file will be rewritten.
output_name: name (input as str) of the csv file for the model results. 
Result figures similar to Figure 6, 7, 8 and 9 in the main text will be plotted at the end of the run (results of monte-carlo simulations are not included). This part can be modified as needed. 
In addition to the key outputs discussed in the main text, the output file also includes calculated sulfur contents at sulfide saturation (SCSS, Smythe et al., 2017), sulfate contents at anhydrite saturation (SCAS, Zajacz and Tsay, 2019) and hydrogen fugacity (fH2, Ohmoto and Kerrick, 1977) along degassing. These values only serve as references if the users are interested, and do not affect the actual degassing calculations. If these results are used, please make sure cite the references listed above accordingly. 
 
1.5 Troubleshooting guides
1.5.1 When testing a new dataset, we suggest the user turn off the fO2 tracker first. This initial run will give a baseline S degassing trend without the impact of evolving fO2 on partition coefficients. 
1.5.2 It is important to check that this initial run has a reasonable S6+/ΣS in the melt for the dataset. If the calculated S6+/ΣS in the melt appears too high or too low, the user should adjust the S-Fe (by changing the parameter S_Fe_choice) relation in the inputs until the model runs with a reasonable initial S6+/ΣS for the given dataset. Then, the user can turn the fO2 tracker on to track the fO2 evolution during degassing. Another way to coarsely decide if the calculated S6+/ΣS is reasonable for the dataset is to compare the measured S concentration to the SCSS and SCAS values. S6+/ΣS could be estimated using SCSS/Smeasured. As we showed in our case studies, if the melt is mostly dominated by S6+, S degassing paths are similar among different S-Fe models while the predicted fO2 varies (e.g. Fuego). If the melt is mostly dominated by S2-, different S6+/ΣS in the melt could result in larger differences in S degassing (e.g. Mauna Kea) and predicted fO2 evolution. 
1.5.3 Large excursions (particularly if the fO2 tracker is enabled) could happen at low pressure degassing due to large partition coefficient and low S in the melt. If Sigma , the tolerance of the calculated log10fO2 (section 1.2.2), is too small or too large, or the pressure step (controlled by parameter l) is too large, it may cause large excursions of fO2 calculations in the run. Therefore, adjusting Sigma and l can help to stabilize the low pressure run. Alternatively, the user can use the monte function to run multiple times and obtain an average results at low pressure degassing. Lastly, if the low-pressure degassing is not of interest (e.g. in studies of comparing to melt inclusion data), the user can simply terminate the degassing run at higher pressure (by change parameter m).
