<p align="center"> 
  <img src="images/fullsizeoutput_2c9.jpeg" alt="Livox-mid-40" width="400px">
</p>
<h1 align="center"> Hydraulic effects of channel reconfiguration and floodplain reconnection </h1>
</br>

<!-- TABLE OF CONTENTS -->
<h2 id="table-of-contents"> :book: Table of Contents</h2>

<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project"> ➤ About The Project</a></li>
    <li><a href="#prerequisites"> ➤ Prerequisites</a></li>
    <li><a href="#Repository Structure"> ➤ Repository Structure</a></li>
    <li><a href="#How to use"> ➤ How to use</a></li>

</details>

![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/aqua.png)

<!-- ABOUT THE PROJECT -->
<h2 id="about-the-project"> :pencil: About The Project</h2>

<p align="justify"> 
This project seeks to develop the evidence base for how channel reconfiguration and floodplain reconnection can affect the hydraulic behavior in the modified reach. This repository provides data and code used to generate figures present within the research article "Hydraulic effects of channel reconfiguration and floodplain reconnection". 
</p>

![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/aqua.png)

<!-- PREREQUISITES -->
<h2 id="prerequisites"> :fork_and_knife: Prerequisites</h2>

**Replicating the outputs** presented in "Hydraulic effects of channel reconfiguration and floodplain reconnection" requires the user to download the data files and code from this GitHub repository, and to be able to run MATLAB 2019a onwards. The easiest way of achieving this is to clone the repository onto you PC. 

![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/aqua.png)

<!-- Repository Structure -->
<h2 id="Repository Structure"> :cactus: Repository Structure</h2>
<p align="justify"> 
  
Below is the an outline of the folder structure within this repository with descriptions provided:
</p>

    .
    ├── code                        # folder containing scripts to reproduce figures
    │   ├── dependencies            # dependencies for running figures scipts
    │   ├── figures                 
    │   │   ├── fig1                # schematic diagram presented in Figure 1
    │   │   ├── fig3                # scripts required to generate outputs presented in Figure 3
    │   │   ├── fig4                # scripts required to generate outputs presented in Figure 4
    │   │   ├── fig5                # scripts required to generate outputs presented in Figure 5
    │   │   ├── fig6                # scripts required to generate outputs presented in Figure 6	
    │   │   ├── fig7                # scripts required to generate outputs presented in Figure 7
    │   │   ├── fig8                # scripts required to generate outputs presented in Figure 8
    │   │   ├── fig9                # scripts required to generate outputs presented in Figure 9
    ├── data                        # folder containing underlying data 
 
  
![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/aqua.png)
  
<!-- How to use -->
<h2 id="How to use"> 👍 How to use</h2>
<p align="justify"> 
    
* Clone or download this repository so that it is accessible on your PC.
* Download the files from: https://data.ncl.ac.uk/articles/dataset/Data_for_output_replication/23501091/2 to your PC.  
* Open MATLAB on your PC.
* To generate Figure 3, run flowCharts.tex from your latex compiling software
* To generate Figure 4, ensure all scripts in "fig4" subfolder are accesible in your MATLAB search path, execute "densityPlot.m", ensuring that you provide the links to the directories containing the relevant datasets (downloaded from Step 2 above).
* To generate Figure 5, ensure all scripts in "fig5" subfolder are accesible in your MATLAB search path, execute "master_fcn_fig5.m", ensuring that you provide the links to the directories containing the relevant datasets (downloaded from Step 2 above).
* To generate Figure 6, ensure all scripts in "fig6" subfolder are accesible in your MATLAB search path, execute "cloud_comparisons.m", ensuring that you provide the links to the directories containing the relevant datasets (downloaded from Step 2 above).
* To generate Figure 7, ensure all scripts in "fig7" subfolders are accesible in your MATLAB search path, execute "cdf_plots.m", ensuring that you provide the links to the directories containing the relevant datasets (downloaded from Step 2 above). 
* To generate Figure 8, ensure all scripts in "fig8" subfolders are accesible in your MATLAB search path, execute "plotGaugingData.m", ensuring that you provide the links to the directories containing the relevant datasets (downloaded from Step 2 above).
* To generate Figure 9, ensure all scripts in "fig9" subfolders are accesible in your MATLAB search path, execute "long_transect_bank_retreat.m", ensuring that you provide the links to the directories containing the relevant datasets (downloaded from Step 2 above).
    
  ![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/aqua.png)

  <p align="center"> 
  <img src="images/IMG_20220128_155022009_HDR.jpg" alt="Goldrill Beck setup" >
  </p>
  <p align="center"> 
  View looking upstream at the Goldrill Beck monitoring site. Shown in the image are the Livox monitoring system, and Riegl VZ-4000 acquiring validation data.
  </p>
  
