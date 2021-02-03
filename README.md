## GRAVEE: Good Risk Assessment Values for Environmental Exposures
### Estimating Point of Departure (POD) from dose-response data using spline meta-regression

GRAVEE is a method for estimating a chemical's Point of Departure (POD) based on
quantitative dose-response data. 

### Input Data  

The input data should be a data frame with two columns, ordered as dose 
and response. Doses should **not** be transformed prior to upload. The
values must be numeric. 

**Example:**  
| Dose  | Response |
| --------- | --------- |
| 1  | 1.4  |
| 1  | 0.28  |
| 1  | 0.9  |
| 5  | 8.32  |
| 5  | 8.6  |
| 5  | 7  |
| 10  | 10.2  |
| 10  | 9.45  |
|...|...|

### Analysis  

All analyses are performed using log<sub>10</sub>-transformed doses.Results
are reported on the original dose scale.  

The core function of this package is `calculate_pod_quantiles()`. GRAVEE 
simulates new samples of dose-response data using bootstrap resampling. Each
simulated sample is fit to a dose-response curve using spline meta-regression.
[Menger Curvature](https://github.com/k-t-to/MengerCurvatureSimulation) is calculated along the interpolated doses of the fitted 
dose-response curve and the POD is reported as the dose with the maximal 
curvature. GRAVEE reports a 95% confidence interval of PODs using the 
calculated PODs from each simulated sample. 
 
### License

This work was created by US government employees working in their government capacity. As a result, this work is in the public domain within the United States. The US Government may exercise copyright in other jurisdictions.

### Contact

Kim T. To - [kimberly.t.to@erdc.dren.mil](kimberly.t.to@erdc.dren.mil)

Lyle Burgoon - [lyle.d.burgoon@erdc.dren.mil](lyle.d.burgoon@usace.army.mil)
