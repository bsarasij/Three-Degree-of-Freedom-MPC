# Three-Degree-of-Freedom-MPC

The 3DoF formulation provides an improved model predictive control (MPC) formulation for. The algorithm relies on a multiple-degree-of-freedom parametrization that enables the user to adjust the speed of setpoint tracking, measured disturbance rejection and unmeasured disturbance rejection independently in the closed-loop system, making it very akin to IMC tuning. Consequently, controller tuning is more flexible and intuitive than relying on objective function weights (such as move suppression) traditionally used in MPC schemes.
<br>
<br>

<p align = "center">
<image src="https://github.com/user-attachments/assets/30a403f3-fbfa-453e-ac14-04cd35d327f8">
</p>
<br>
<br>

The 3DoF is obtained by a set of filters, the first one being for setpoint tracking. This is typically a type-I discrete-time filter with a tuning parameter $\alpha_r$ that adjusts the output response for setpoint tracking. The second degree of freedom arises from the measured disturbance filter whose speed of response $\alpha_d$ dictates how aggressively the controller reacts to such dist, This filter is typically type-I as well but can also be type-II depending on the system requirements. 
<br>
<br>
<p align="center">
<img src="https://github.com/user-attachments/assets/e2154192-d6df-4abd-8225-a27770ea5e88">
</p>
<br>
<br>

For the 3rd degree we use a set of two Kalman filters that perform state estimation while decoupling the effects of measured and unmeasured disturbances on the output.
<p align="center">
<img src="https://github.com/user-attachments/assets/8d318034-0fe3-4de6-bbe4-b94f9d2b3f96">
</p>

<br>
<br>

## Application demonstrated using a microalgae photobioreactor

<br>
There is a semi-industrial raceway photobioreactor facility located at the IFAPA Research Center near the UAL campus in Andalusia, Spain. They perform combined wastewater treatment and microalgae biomass production for addressing the key challgenge of greenhouse gas management in a sustainable manner. The idea is to regulate the CO2 inflow into the reactor to maintain the pH at the optimal value of 8 in presence of changing measured and unmeasured disturbances such as Temperature (measured), Radiation (measured) , Airflow (unmeasured), Dilution (unmeasured).
<br>
<br>
<p align="center">
<img src="https://github.com/user-attachments/assets/1379b1b3-4a89-45f2-b2a5-9b0914d698fc">
</p>

<br>
When applied to the photobioreactor system for a constant setpoint of pH = 8, the 3DoF controller demonstrates: 

* Setpoint tracking of pH during the initial part of the day with no offset or oscillation. 

* Active rejection of measured Radiation and Temperature disturbances which change throughout the day. 

* Rapid mitigation of noise and robustness to plant-model mismatch. The controller also successfully rejects unmeasured disturbances arising from Dilution and airflow changes.

* Physically realizable CO2 inflow into the system that demonstrates plant-friendliness.  

<br>
<br>
<p align="center">
<img src="https://github.com/user-attachments/assets/2c4726c2-7c65-46e5-bb19-b8bce2dd94e2">
</p>




