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
<img src="https://github.com/user-attachments/assets/22da636d-8d91-478c-b5a1-d3f46e56383a">
</p>
