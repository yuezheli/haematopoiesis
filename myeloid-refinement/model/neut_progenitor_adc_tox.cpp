[PROB]

To fit for ADC tox on neutrophil prognitor from Zhao et al., 2017
https://aacrjournals.org/mct/article/16/9/1866/146539/A-Potential-Mechanism-for-ADC-Induced-Neutropenia
  
Neutrophil development model structure & parameters obtained from https://pubmed.ncbi.nlm.nih.gov/28303306/ & https://immunenetwork.org/DOIx.php?id=10.4110/in.2017.17.5.298

[SET]
delta = 0.1, end = 6

[CMT] 
promyelocyte
myelocyte01
myelocyte02
metamyelocyte
bandcell
neutrophil

cyto_ADC

[PARAM]

betaGM = 2.5 // GM death rate; unit day-1; estimated from half life of 6.7h; https://pubmed.ncbi.nlm.nih.gov/28303306/
a_myelocyte = 32  // back calculate from myelocyte -> metamyelocyte transition time; https://www.nature.com/articles/nature14242
tau_promyelocyte = 2.08 // mean residence time of proliferating neutrophil progenitors, [day], https://pubmed.ncbi.nlm.nih.gov/28303306/
tau_myelocyte = 2.08 // mean residence time of proliferating neutrophil progenitors, [day], https://pubmed.ncbi.nlm.nih.gov/28303306/
tau_metamyelocyte = 0.5 // mean residence time for nonproliferating progenitor cells for neutrophil (metamyelocyte & band cell), [day],

medium_ADC = 0 // assumed to be a constant

k_in_ADC = 3 // tuned 

n_ADC_neut = 1.13 // parmaeters from T-DM1 on CFU-MK cells

// IC50 for ADC-induced neutrophil progenitor death, [nM], https://aacrjournals.org/mct/article/16/9/1866/146539/A-Potential-Mechanism-for-ADC-Induced-Neutropenia
IC50_ADC = 14  


[MAIN]
double amp_myelocyte = a_myelocyte * pow(IC50_ADC, n_ADC_neut) / ( pow(IC50_ADC, n_ADC_neut) + pow(cyto_ADC, n_ADC_neut) ); 
  
double tau1_myelocyte = (log2(amp_myelocyte)-1)/log2(amp_myelocyte) * tau_myelocyte;
double tau2_myelocyte = tau_myelocyte - tau1_myelocyte;
double krep1_myelocyte = (amp_myelocyte/2-1)/ tau1_myelocyte;
double krep2_myelocyte = 1/ tau2_myelocyte;
double k12_myelocyte = amp_myelocyte/(2*tau1_myelocyte);
double k_myelocyte2metamyelocyte = 2/tau2_myelocyte;

[ODE]
dxdt_cyto_ADC = k_in_ADC * (medium_ADC - cyto_ADC) - k12_myelocyte * cyto_ADC; 

dxdt_promyelocyte = - 1/tau_promyelocyte * promyelocyte; 

dxdt_myelocyte01 = 2 * 1/tau_promyelocyte * promyelocyte + (krep1_myelocyte - k12_myelocyte) * myelocyte01; 
dxdt_myelocyte02 = k12_myelocyte * myelocyte01 + (krep2_myelocyte - k_myelocyte2metamyelocyte) * myelocyte02; 

dxdt_metamyelocyte = k_myelocyte2metamyelocyte * myelocyte02 - 1/tau_metamyelocyte * metamyelocyte;
dxdt_bandcell = 1/tau_metamyelocyte * metamyelocyte - 1/tau_metamyelocyte * bandcell;
dxdt_neutrophil = 1/tau_metamyelocyte * bandcell - betaGM * neutrophil;

[TABLE]
capture totalcell = promyelocyte + myelocyte01 + myelocyte02 + metamyelocyte + bandcell + neutrophil;

[CAPTURE]
amp_myelocyte, k12_myelocyte