# Mathematical Models in Asteroid Impact Risk Assessment

[cite_start]This document outlines the key mathematical equations used in the Probabilistic Asteroid Impact Risk (PAIR) model to assess the threat from sub-300 m asteroids[cite: 6, 28]. [cite_start]The model incorporates physics-based equations for atmospheric entry, breakup, and ground damage estimation within a Monte Carlo framework[cite: 29, 88].

***

## Atmospheric Entry and Breakup Modeling ☄️

[cite_start]The model simulates an asteroid's flight through the atmosphere using a set of single-body ordinary differential equations to track its changing mass, velocity, and trajectory[cite: 109].

### Flight Dynamics Equations

[cite_start]These equations are integrated starting from an altitude of 100 km[cite: 129]:

$$\frac{dm}{dt} = -0.5 \rho_{air} v^3 A \sigma \quad \text{(1a)}$$
$$\frac{dv}{dt} = -0.5 \rho_{air} v^2 A C_D / m - g \sin\theta \quad \text{(1b)}$$
$$\frac{d\theta}{dt} = \left(\frac{v}{R_E+h} - \frac{g}{v}\right) \cos\theta \quad \text{(1c)}$$
$$\frac{dh}{dt} = v \sin\theta \quad \text{(1d)}$$

* [cite_start]**Variables Defined**[cite: 128]:
    * **$m$**: Mass of the asteroid.
    * **$v$**: Velocity of the asteroid relative to the atmosphere.
    * **$\theta$**: Flight path angle.
    * **$h$**: Altitude.
    * **$t$**: Time from the initial atmospheric interface.
    * **$g$**: Acceleration due to gravity.
    * **$\rho_{air}$**: Local atmospheric density.
    * **$R_E$**: Radius of the Earth.
    * **$A$**: Cross-sectional area of the object.
    * **$C_D$**: Drag coefficient.
    * **$\sigma$**: Ablation coefficient.

---

### Fragmentation and Dispersion

[cite_start]Fragmentation begins when the aerodynamic pressure on the asteroid exceeds its strength[cite: 129].

[cite_start]**Breakup Condition**[cite: 130]:
$$\rho_{air} v^2 > \text{aerodynamic strength} \quad \text{(2)}$$

[cite_start]**Fragment Strength Scaling**[cite: 134]:
After the initial breakup, the strength of the resulting "child" fragments increases based on a power law.
$$S_{child} = S_{parent} \left( \frac{m_{parent}}{m_{child}} \right)^a \quad \text{(3)}$$

* [cite_start]**Variables Defined**[cite: 136]:
    * **$S$**: The aerodynamic strength of the parent or child fragment.
    * **$m$**: The mass of the parent or child fragment.
    * **$a$**: The power law strength-scaling exponent.

[cite_start]**Debris Cloud Dispersion**[cite: 140, 141]:
The cloud of debris created during fragmentation spreads out as it travels through the atmosphere.
$$v_{dispersion} = v_{cloud} \left( \frac{3.5 \rho_{air}}{\rho_{cloud}} \right)^{1/2} \quad \text{(4)}$$

* [cite_start]**Variables Defined**[cite: 143]:
    * **$v_{cloud}$**: The velocity of the debris cloud.
    * **$\rho_{cloud}$**: The density of the initial asteroid material.
    * **$\rho_{air}$**: Local atmospheric density.

***

## Damage Modeling 💥

[cite_start]The model estimates ground damage from two primary hazards: blast overpressure and thermal radiation[cite: 202]. [cite_start]The larger of the two damage areas is used for the final assessment[cite: 202].

### Blast Overpressure Damage

[cite_start]The radius of ground damage is calculated using a scaling relation based on nuclear test data, with a baseline threshold of **4 psi overpressure**[cite: 163, 172].

$$R_{ground} = 2.09h - 0.449h^2 E^{-1/3} + 5.08E^{1/3} \quad \text{(5)}$$

* [cite_start]**Variables Defined**[cite: 179]:
    * **$R_{ground}$**: The ground damage radius in kilometers.
    * **$h$**: The burst altitude in kilometers.
    * **$E$**: The impact energy in megatons (Mt).

---

### Thermal Radiation Damage

[cite_start]The model estimates the radius where thermal radiation is intense enough to cause **3rd-degree burns**[cite: 181, 191].

[cite_start]**Threshold Radius**[cite: 188]:
$$r = \sqrt{\frac{\eta E}{2\pi \Phi_j}} \quad \text{(6)}$$

* [cite_start]**Variables Defined**[cite: 189, 190]:
    * **$r$**: The threshold radius from the burst point.
    * **$\eta$**: Luminous efficiency, the fraction of energy released as thermal radiation.
    * **$E$**: The impact energy.
    * **$\Phi_j$**: The thermal exposure (heat per unit area) for a specific damage level.

[cite_start]**Ground Radius for Airbursts**[cite: 196]:
For an airburst, this spherical radius is projected onto the ground to find the affected area.
$$R_{ground} = \sqrt{r^2 - h^2}$$

* [cite_start]**Variables Defined**[cite: 196]:
    * **$h$**: The burst altitude.

***

## Input Parameter Calculations 🎲

[cite_start]These equations are used within the model's Monte Carlo framework to generate a diverse set of impact scenarios by sampling from statistical distributions[cite: 90].

### Asteroid Diameter

[cite_start]The asteroid's diameter is calculated from its observable properties: **absolute magnitude ($H$)** and **albedo ($p_v$)**[cite: 244, 246].

$$D = (1.326 \times 10^6) \times 10^{-H/5} / \sqrt{p_v}$$

* [cite_start]**Variables Defined**[cite: 247]:
    * **$D$**: Diameter in meters.
    * **$H$**: The absolute magnitude.
    * **$p_v$**: The albedo.

---

### Entry Angle

[cite_start]The entry angle is sampled from a distribution where a **45° entry is most probable**, and vertical or very shallow entries are much less likely[cite: 270].

$$\theta^\circ = (90^\circ/\pi) \cos^{-1}(2U-1)$$

* [cite_start]**Variables Defined**[cite: 271]:
    * **$\theta^\circ$**: The entry angle in degrees.
    * **$U$**: A random number sampled uniformly between 0 and 1.
