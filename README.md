This project aims to reproduce the codes from the paper by Schmalholz et al. (2024), using the Julia language.
Below, I briefly describe the different programs. For more details, you can refer to the comments within the code or the dedicated documentation section.

1 – HM_model_I.jl
This program solves the diffusion equations for fluid pressure ($Pf$) and porosity ($\phi$).
The equations have the same structure for both variables and include poroelastic effects (Yarushina & Podladchikov, 2015).
(De)hydration reactions are not considered in this model.
This simplified case is used to test and compare various numerical methods that I am familiar with, including:
- Explicit method
- Implicit method
- Pseudo-transient method
