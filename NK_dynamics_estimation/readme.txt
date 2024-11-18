NOTE NK CELL DYNAMICS IS NOT INCLUDED IN THE FINAL MODEL. 
THIS FOLDER IS USED AS A PLACEHOLDER FOR POTENTIAL NK CELL DYNAMICS DEVELOPMENT.

Goal: 

The goal of this folder is to estimate the dynamics of NK cells in mouse.


Source of data: 
https://mcb.berkeley.edu/labs/raulet/Resources/04-01%20NKDynamicsJI.pdf


Equation: 
dNK/dt = s + alpha * NK - beta * NK
where s is the differentiation from progenitor cells (not explicitly modeled here), alpha is the proliferation rate (unit: day-1) and beta is the death rate (unit: day-1). 

By solving this ODE, then we get NK = C0 + s * t + exp( (alpha - beta) * t), where C0 is a constant. 

Fitting: 

Scenario in Fig 2A, i.e. NK cells are continuously labeled. This gives the differentiation rate from progenitor cells. The exponential part (i.e. NK cell division and death) cannot be captured in the fitting. 
    s = 43655, i.e. NKP * differentiation_rate = 43655

Scenario in Fig 3A, i.e. NK cell death rate in spleen; rate = 0.025 day-1

Scenario ib Fig 3B, i.e. NK cell death rate in bone marrow; rate = 0.0072 day-1