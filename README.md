# Scramjet Isolator

This project focuses on the **analytical and computational study of scramjet isolator flows**, where the effect of changing the **equivalence ratio (ϕ)** on flow properties is investigated. The study is based on governing equations from *RTO-EN-AVT-185 (Scramjet Isolators – Michael K. Smart)* and solved numerically in **MATLAB** using the `ode45` solver.  

## Project Overview
The scramjet isolator plays a critical role in dual-mode scramjet operation, ensuring stable combustion and preventing inlet unstart. In this project:  
- A **1D analytical model** was developed to capture variations in **Mach number, pressure, temperature, and area ratios** along the isolator–combustor duct.  
- Governing equations account for **friction, heat addition, and area change** in the flow.  
- The **equivalence ratio (ϕ)** was varied between **0 and 1**, and its impact on flow properties was analyzed.  
- Results were visualized in the form of **graphs** showing trends of Mach number, pressure ratio (P/P₂), temperature ratio (T/T₂), and area ratios (A/A₂, Ac/A₂).  

## Key Results
- **ϕ = 0.5** → Attached supersonic flow with stable isolator operation.  
- **ϕ = 0.72** → Shock-induced separation and **supersonic reattachment** observed.  
- **ϕ = 0.81** → Formation of a **thermal throat**, leading to subsonic reattachment in the isolator.  

These results align with published trends, demonstrating the model’s ability to capture **dual-mode scramjet behavior**.  
