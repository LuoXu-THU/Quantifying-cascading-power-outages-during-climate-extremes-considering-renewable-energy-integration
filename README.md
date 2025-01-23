# Code for "Quantifying cascading power outages during climate extremes considering renewable energy integration"

The code is developed for paper 'Quantifying cascading power outages during climate extremes considering renewable energy integration'

This repository contains the MATLAB code developed for simulating spatiotemporal power outages under climate extremes, specifically focused on a case study of the Puerto Rico power grid.


## Installation Requirements

The code is running in MATLAB R2022a. Install Gurobi 10.1, YALMIP R20230609, MATPOWER 7.0.


## Functionality

The methodology for simulating the wind field of Hurricane Fiona in 2022 draws upon the foundational work of Chavas, Lin, and Emanuel (2015), who developed a comprehensive model to depict the full radial structure of the tropical cyclone wind field. This model provides a physics-based method to simulate hurricane wind profiles. The reference for this method is as follows:

\* Chavas, D., N. Lin, and K. Emanuel (2015). A model for the complete radial structure of the tropical cyclone wind field. Part I: Comparison with observed structure. J. Atmos. Sci. 

The wind field simulation for Hurricane Fiona utilizes input data from the International Best Track Archive for Climate Stewardship (IBTrACS) database, which offers comprehensive global tropical cyclone tracking information. See:

\* Schreck, C. J., Knapp, K. R., & Kossin, J. P. (2014). The impact of best track discrepancies on global tropical cyclone climatologies using IBTrACS. *Monthly Weather Review*, 142, 3881-3899.

'Cascading_power_outages.m' serves as the main function to generate spatiotemporal power outages.



## Contact Information

For further information or inquiries, please contact Luo Xu (luoxu@princeton.edu)



## License

Copyright © 2025, Luo Xu. All rights reserved. When the paper is published, released under the Creative Commons Attribution 4.0 International License (CC BY 4.0). 