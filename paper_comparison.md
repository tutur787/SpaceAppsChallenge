I am comparing whitepapers with my team for this hackathon:
- Evaluate which paper makes more sense
- See if we can grab different resources from the whitepapers to cover all aspects we need
We need kind of like a blueprint of the math:
<hackathon_description>
# 2025 NASA Space Apps Challenge
## Meteor Madness
Summary
A newly identified near-Earth asteroid, "Impactor-2025," poses a potential threat to Earth, but do we have the tools to enable the public and decision makers to understand and mitigate its risks? NASA datasets include information about known asteroids and the United States Geological Survey provides critical information that could enable modeling the effects of asteroid impacts, but this data needs to be integrated to enable effective visualization and decision making. Your challenge is to develop an interactive visualization and simulation tool that uses real data to help users model asteroid impact scenarios, predict consequences, and evaluate potential mitigation strategies.

Background
The discovery of near-Earth asteroids like "Impactor-2025" highlights the ongoing risk of celestial objects colliding with our planet, potentially causing catastrophic damage. Asteroid impacts, though rare, could cause widespread devastation, including tsunamis, seismic events, and atmospheric changes. NASA’s Near-Earth Object (NEO) program tracks thousands of asteroids and data from NASA’s NEO Application Programming Interface (API) provides asteroid characteristics (e.g., size, velocity, orbit). The U.S. Geological Survey (USGS) offers environmental and geological datasets (e.g., topography, seismic activity, tsunami zones) critical for modeling impact effects. However, these datasets are often siloed; the ability to predict and mitigate specific impact scenarios remains limited by the complexity of integrating these diverse datasets and communicating risks to stakeholders.

Existing tools lack user-friendly interfaces to simulate and visualize impacts for scientists, policymakers, and the public. In addition, current tools often focus on detection and orbital tracking but fall short in simulating impact consequences or evaluating mitigation strategies like deflection (e.g., kinetic impactors or gravity tractors). These tools are also often either too technical for public use or overly simplistic, missing key environmental impacts. Such gaps hinder preparedness efforts and public understanding of potential impacts if an asteroid were to collide with Earth. A tool that combines accurate data integration, realistic physics-based simulations, and intuitive visualizations could bridge the gap between complex science and actionable insights.

Objectives
Your challenge is to develop an interactive visualization and simulation tool that enables users to explore asteroid impact scenarios, predict consequences, and evaluate mitigation strategies using real NASA and USGS datasets. Can you create a web-based platform that integrates asteroid parameter data (e.g., size, velocity, trajectory) from NASA’s NEO API with USGS datasets pertaining to environmental and geological impacts (e.g., tsunami zones, seismic activity, topography), and transforms that raw data into a powerful educational and decision-support tool for global asteroid risk management?

Think about how your tool can incorporate intuitive controls and dynamic visualizations, such as animated orbital paths or impact zone maps, to engage scientists, policymakers, and the public. How will your tool simulate the asteroid’s trajectory using orbital mechanics, calculate energy from an impact (e.g., determine the crater size, seismic magnitude), and provide visualizations of the outcomes of mitigation strategies? Your tool could allow users to input variables (e.g., asteroid size and velocity) or specify deflection strategies (e.g., kinetic impactors) and provide visualizations of the outcomes, including the asteroid’s path relative to Earth, potential impact points, and environmental effects like tsunamis or atmospheric changes.

Consider how your tool can prioritize scientific accuracy, user-friendliness, and educational value. For example, you could potentially incorporate explanatory tooltips or infographics into your tool to demystify complex concepts. Be creative! Can your tool incorporate gamified or storytelling elements to engage users, or leverage machine learning to predict impact outcomes—all while ensuring scientific accuracy and accessibility? For instance, you could design a scenario where users “defend” Earth by adjusting deflection parameters.

How will your tool balance scientific rigor with accessibility, making complex astronomical threats understandable and actionable for all? Feel free to explore use of various formats, such as a web application with 3D visualizations of asteroid trajectories and impact zones, an interactive dashboard that allows the user to adjust variables like asteroid size or deflection timing, or a gamified simulation. But don’t forget to focus on leveraging NASA and USGS data to tell a compelling story about asteroid threats. For example, your tool could simulate a “what-if” scenario for Impactor-2025, showing how a small change in velocity alters the impact point, or map tsunami risks using USGS elevation data.

This challenge empowers you to transform raw data into a powerful educational and decision-support tool for global asteroid risk management!

Potential Considerations
You may (but are not required to) consider the following:

General Guidance

Target Audience: Scientists, policymakers, educators, and the public could use your tool, so try to ensure accessibility for non-experts while retaining technical depth.
Scalability: You could build a modular system to handle additional datasets (e.g., atmospheric density, population density).
Performance: How can you optimize simulations and visualizations for smooth browser performance, especially for 3D rendering?
Execution: You are encouraged to use technologies for backend data processing and for user interfaces.
Scientific Considerations

Orbital Mechanics: You can model the asteroid’s trajectory using Keplerian orbital elements (e.g., semi-major axis, eccentricity, true anomaly) with standard orbital position calculations.
Impact Energy: You can estimate kinetic energy based on the asteroid’s mass (derived from size and density, e.g., 3000 kg/m³) and velocity, then convert to the Trinitrotoluene (TNT) equivalent for impact scale.
Crater Scaling: You can use established scaling relationships to estimate crater size based on impact energy.
Environmental Effects: You can leverage USGS data to model secondary effects like tsunamis (using coastal elevation) or seismic waves (for inland impacts).
Technical Tips

Technologies: You could use Python (Flask/Django) for backend processing, JavaScript (Three.js/D3.js) for visualizations, and HTML/CSS for interfaces.
Interactivity: Your tool could include sliders, dropdowns, or maps for user inputs with real-time visualization updates.
Visualization: You could use 3D for orbital paths (Three.js) and 2D for impact maps (D3.js).
Error Handling: Don’t forget to implement fallbacks for potential API failures during the hackathon.
Pitfalls to Avoid

Overcomplication: Consider avoiding use of complex physics models (e.g., n-body simulations) that slow the tool or exceed hackathon constraints.
Data Misuse: Take care to correctly interpret NASA and USGS data and verify units (e.g., NEO API’s miss distance is in kilometers).
Non-Intuitive User Interface: Try to avoid cluttered interfaces or technical jargon; don’t forget to test for user clarity.
Ignoring Mitigation: You could include deflection strategies (e.g., velocity changes via kinetic impactors) to show proactive solutions.
Standout Features

Gamification: You could create a “defend Earth” mode where users test deflection strategies under time pressure.
Educational Overlays: You might add tooltips or pop-ups explaining terms like “eccentricity” or “impact energy.”
Regional Focus: Consider allowing the user to zoom into specific regions (e.g., coastal cities) for localized impact visualizations.
Mitigation Scenarios: You could simulate advanced deflection methods (e.g., gravity tractors, laser ablation).
Storytelling: You could frame the tool as an interactive narrative, guiding users through a hypothetical Impactor-2025 scenario.
Accessibility: Can you include colorblind-friendly palettes, keyboard navigation, and multilingual support?
Add-Ons to Consider

Real-Time Data: Your tool could fetch live NEO data from NASA’s API.
Social Sharing: You could allow users to share simulation results (e.g., impact maps) on social media.
Mobile Compatibility: Will you optimize the tool for use on mobile browsers?
Augmented Reality (AR): You could explore using AR frameworks (e.g., A-Frame) to project asteroid paths in real-world environments.
For data and resources related to this challenge, refer to the Resources tab at the top of the page.

https://api.nasa.gov/
NASA Near-Earth Object (NEO) Web Service Application Programming Interface (API): This API provides access to NASA’s Near-Earth Object (NEO) dataset, including orbital parameters (e.g., semi-major axis, eccentricity, inclination), size estimates, and close-approach data for asteroids and comets. Participants can use this data to model the trajectory of a hypothetical asteroid like Impactor-2025, calculate potential impact probabilities, and simulate deflection strategies. The API is publicly accessible with a free key (directions to generate a key are on the above website), making it ideal for real-time data integration into web-based tools. It’s particularly useful for generating realistic orbital paths and impact scenarios. After visiting the above API link, scroll down and open the “Asteroids - NeoWs” section for more information on how to use the API.


U.S. Geological Survey (USGS) National Earthquake Information Center (NEIC) Earthquake Catalog: The USGS NEIC provides a comprehensive, publicly accessible catalog of global earthquake data, including location, magnitude, and depth. Participants can use this dataset to model seismic effects of an asteroid impact by correlating impact energy with equivalent earthquake magnitudes. The data can also inform visualizations of ground shaking or secondary effects in specific regions, enhancing the ability to depict realistic environmental consequences. The dataset’s open access and global coverage make it valuable for international teams. Learn more about USGS Catalog Event Terms here.
https://earthquake.usgs.gov/earthquakes/search/

USGS National Map Elevation Data: This dataset offers high-resolution digital elevation models (DEMs) for the United States and some global regions, available through the USGS Earth Resources Observation and Science (EROS) Center. Participants can use elevation data to model tsunami inundation or crater formation in specific impact locations, enhancing visualizations of environmental impacts. The data is publicly available, downloadable in formats like GeoTIFF, and suitable for integration into Geographic Information System (GIS)-based or 3D visualization tools, making it essential for simulating localized geological effects of an asteroid strike.

https://www.usgs.gov/programs/national-geospatial-program/national-map

USGS National Map Training Videos: The National Map Training page from the U.S. Geological Survey offers a rich collection of instructional videos designed to help users understand and utilize The National Map (TNM) products and services. These resources are tailored for geospatial professionals, scientists, and the general public.
https://www.usgs.gov/programs/national-geospatial-program/national-map-training


Approximate Positions of the Planets: This NASA resource provides Keplerian orbital parameters for planets in our solar system and formulas for calculating their positions. While focused on planets, the Keplerian parameters (e.g., semi-major axis, eccentricity) are applicable to asteroid orbit modeling. Participants can use this resource as a reference to understand orbital mechanics and adapt the principles to simulate the trajectory of Impactor-2025 or other asteroids, making it a foundational resource for accurate orbital visualizations.
https://ssd.jpl.nasa.gov/planets/approx_pos.html

Small-Body Database Query Tool: NASA’s Small-Body Database Query tool provides access to Keplerian parameters for near-Earth objects (NEOs), including potentially hazardous asteroids (PHAs). Participants can query specific asteroid data (e.g., size, orbit, velocity) to model the trajectory and impact scenarios for Impactor-25 and other asteroids. This publicly accessible database is critical for obtaining real asteroid data, enabling realistic simulations and visualizations in the challenge. Other Small-Body Database Tools may also be useful
https://ssd.jpl.nasa.gov/tools/sbdb_query.html

Near-Earth Comets - Orbital Elements API: This NASA Open Data Portal API provides Keplerian orbital elements for 170 near-Earth comets (NECs) in JavaScript Object Notation (JSON) format. While focused on comets, the orbital parameters are similar to those for asteroids, making it a useful resource for participants to practice retrieving and processing orbital data. It supports the development of orbital propagators for Impactor-2025, enhancing the ability to simulate trajectories.
https://data.nasa.gov/dataset/near-earth-comets-orbital-elements-api

Elliptical Orbit Simulator: This NASA tutorial explains how to simulate and plot elliptical orbits using Keplerian parameters, with examples in R. Participants with intermediate coding skills can reuse or port the orbital propagator code to model the paths of Impactor-2025 or other asteroids. The tutorial is publicly accessible and provides a clear starting point for trajectory simulations, with a reminder to cite the source if used.
https://nasa.github.io/mission-viz/RMarkdown/Elliptical_Orbit_Design.html

Eyes on Asteroids - NASA/Jet Propulsion Laboratory (JPL): This NASA/JPL static orrery visualizes orbital trajectories and labeled asteroids in 3D. Participants can study its design to create similar visualizations for Impactor-2025 and other asteroids, focusing on user-friendly displays of asteroid paths and impact risks. The public accessibility and focus on asteroids make it a valuable reference for interactive visualization ideas.
https://eyes.nasa.gov/apps/asteroids/

Near-Earth Object Surveillance Satellite (NEOSSAT) - Astronomy Data: The dataset includes the astronomical images from the Near-Earth Object Surveillance Satellite (NEOSSat). NEOSSat is the world's first space telescope dedicated to detecting and tracking asteroids, comets, satellites and space debris.
https://donnees-data.asc-csa.gc.ca/en/dataset/9ae3e718-8b6d-40b7-8aa4-858f00e84b30

Near-Earth Object Surveillance Satellite (NEOSSAT): observing asteroids, space debris and exoplanets: This website provides more information on the Near - Earth Object Surveillance Satellite (NEOSSat). It covers topics such as asteroid and space debris tracking, exoplanet detection, and Canada’s contributions to space situational awareness.
https://www.asc-csa.gc.ca/eng/satellites/neossat/

</hackathon_description>

<A probabilistic asteroid impact risk model: assessment of sub-300 m
impacts>
Icarus 289 (2017) 106–119
Contents lists available at ScienceDirect
Icarus
journal homepage: www.elsevier.com/locate/icarus
A probabilistic asteroid impact risk model: assessment of sub-300 m
impacts
Donovan L. Mathias a,∗
, Lorien F. Wheeler b
, Jessie L. Dotsonc
a NASA Ames Research Center, MS 258-5, Moffett Field, CA 94035, United States b CSRA, NASA Ames Research Center, MS 258-6, Moffett Field, CA 94035, United States c NASA Ames Research Center, MS 244-30 Moffett Field, CA 94035, United States
a r t i c l e i n f o
Article history:
Received 1 November 2016
Revised 24 January 2017
Accepted 15 February 2017
Available online 20 February 2017
Keywords:
Asteroids
Near-Earth objects
Impact processes
a b s t r a c t
A comprehensive asteroid threat assessment requires the quantification of both the impact likelihood
and resulting consequence across the range of possible events. This paper presents a probabilistic asteroid impact risk (PAIR) assessment model developed for this purpose. The model incorporates published
impact frequency rates with state-of-the-art consequence assessment tools, applied within a Monte Carlo
framework that generates sets of impact scenarios from uncertain input parameter distributions. Explicit
treatment of atmospheric entry is included to produce energy deposition rates that account for the effects of thermal ablation and object fragmentation. These energy deposition rates are used to model the
resulting ground damage, and affected populations are computed for the sampled impact locations. The
results for each scenario are aggregated into a distribution of potential outcomes that reflect the range of
uncertain impact parameters, population densities, and strike probabilities. As an illustration of the utility
of the PAIR model, the results are used to address the question of what minimum size asteroid constitutes a threat to the population. To answer this question, complete distributions of results are combined
with a hypothetical risk tolerance posture to provide the minimum size, given sets of initial assumptions
for objects up to 300 m in diameter. Model outputs demonstrate how such questions can be answered
and provide a means for interpreting the effect that input assumptions and uncertainty can have on final
risk-based decisions. Model results can be used to prioritize investments to gain knowledge in critical
areas or, conversely, to identify areas where additional data have little effect on the metrics of interest.
Published by Elsevier Inc.
This is an open access article under the CC BY-NC-ND license.
(http://creativecommons.org/licenses/by-nc-nd/4.0/)
1. Introduction
Dramatic evidence suggests extreme consequences can result
when asteroids strike the earth. Images of tree-fall in Tunguska
(Vasilyev, 1998), the impressive Meteor Crater (Barringer, 1909),
historical record of the Eltanin ocean impact (Gersonde et al.,
1997), and the K-T dinosaur extinction event (Alvarez et al., 1980)
all represent potential impact consequences. In the 1990s, studies began to quantify the level of hazard posed by such asteroid
strikes by evaluating their potential consequences and expected
impact rate (Chapman and Morrison, 1994; Morrison et al., 1994;
Toon et al., 1997). In 1992 the threat was considered significant
enough that the U.S. Congress mandated the Spaceguard Survey to
locate 90% of near earth asteroids greater than 1 km in diameter
(Morrison, 1992). The next major risk assessment was produced in
∗ Corresponding author.
E-mail addresses: donovan.mathias@nasa.gov (D.L. Mathias),
lorien.wheeler@nasa.gov (L.F. Wheeler), jessie.dotson@nasa.gov (J.L. Dotson).
2003 (Stokes et al., 2003) and in 2005 led to the George E. Brown
Jr. survey to extend the search criteria down to objects 140 m in
diameter (NASA Authorization Act, 2005). In 2010, the National Research Council (NRC) produced a report (National Research Council, 2010), which asserted that a search for much smaller objects
was required. Then in 2013, while the world watched the flyby
of asteroid 2012 DA14, a smaller asteroid only 20 m in diameter—
below the threshold many considered to be a threat—entered the
atmosphere over Chelyabinsk, Russia and caught the attention of
the world with widespread video capture of the resulting fireball.
While the Chelyabinsk meteor caused no fatalities, the shockwave
produced by its atmospheric breakup injured over 1000 people and
caused $33 million in damage (Popova et al., 2013). Following the
Chelyabinsk event, a number of new looks at impact risk assessment have emerged (Rumpf et al., 2016; Motiwala et al., 2015). In
the current paper, further advances to impact risk modeling are
presented.
When reviewing the literature, it is worth differentiating between the impact hazard and the impact risk. A hazard is typhttp://dx.doi.org/10.1016/j.icarus.2017.02.009
0019-1035/Published by Elsevier Inc. This is an open access article under the CC BY-NC-ND license. (http://creativecommons.org/licenses/by-nc-nd/4.0/)
D.L. Mathias et al. / Icarus 289 (2017) 106–119 107
ically considered to be a potential damage-causing event or action (Ericson, 2005), whereas risk is defined as the consequence
of an event weighted by the probability of the event occurring.
By this convention, airburst blast damage would be considered a
class of hazard, while the expected casualties (average number of
casualties per year) would be a measure of the risk because it includes both the consequence of the impacts as well as the likelihood of such events. There are different interpretations among
other sources, but we have adopted this convention to remain consistent with the risk assessment community and with recent papers such as Reinhardt et al. (2016) and Rumpf et al. (2016).
Previous papers have looked at impact consequences for representative events (Melosh, 2007; Boslough and Crawford, 2007;
Ivanov et al., 1997; Kring, 1997; Chyba et al., 1993) and collectively
from a hazard perspective (Toon et al., 1997; Chapman and Morrison, 1994; Morrison et al., 1994). Risk assessments have either
evaluated the risk posed by known objects (Chesley et al., 2002;
Rumpf et al., 2016) or by a statistical ensemble of hypothetical objects (Stokes et al., 2003; National Research Council, 2010; Chapman, 2004; Reinhardt et al., 2016; Motiwala et al., 2015). Stokes et
al. performed the first real look at quantifying risk by considering
a range of impactor sizes, each with a single, representative set of
physical properties (density, velocity, strength, etc.). Subsequent efforts have extended Stokes’s approach (National Research Council,
2010; Chapman, 2004), but the updates largely focus on improvements in characterization of the asteroid population and corresponding impact frequency estimates. The Stokes et al. (2003), National Research Council (2010), and Chapman (2004) assessments
all considered representative impact scenarios with an assumed
set of mean asteroid and entry parameters (i.e., density, entry angle, velocity, etc.) and, often, average human population distribution on Earth. As a result, the metrics used to report the risk have
typically been based on annual expected casualties, i.e., the average annual casualty rate calculated based on infrequent events averaged over a long time period. Stokes et al. presented expected
casualty rates, but also considered the actual population distribution and presented probability numbers describing the likelihood
of affecting more than a given number of people as a function of
impactor size, assuming a single representative scenario. Reinhardt
et al. (2016) and Motiwala et al. (2015) both included uncertainty
distributions that described the range of possible impactor characteristics as well as random impact locations. These authors both
performed the risk assessment using a Monte Carlo approach and
presented uncertainty distributions of outcomes. Motiwala et al.
(2015) included the atmospheric entry and breakup process explicitly in the risk modeling, whereas Reinhardt utilized parametric
equations to estimate the ground hazards for each scenario.
The current work extends the approach implemented in
Motiwala et al. (2015) and embeds physical meteor entry models into a probabilistic risk assessment. Single-body equations are
explicitly integrated for each impact case within the Monte Carlo
framework so that a more faithful treatment of entry physics can
be associated with the consequence assessment. This allows explicit representation of the physical processes that occur during
entry, increasing the fidelity of the results, and allows for the examination of the sensitivity of the results to the input assumptions.
The following sections will describe the risk model framework and
the physical models currently incorporated. While the model can
be used to assess the threat of impact in a variety of scenarios,
including the case of a specific object, this paper shows ensemble risk results for a sample range of stochastic impact scenarios.
These results include the distributions of damage likelihoods given
the input assumptions, and illustrate how the distributions provide much more information than average expected values. As an
example, the results are used to define a minimum asteroid size
that would constitute a considerable threat, based on a hypothetical risk tolerance. Finally, we demonstrate how the current risk
model can be used to examine how this hypothetical threat size
changes based on the assumptions and uncertainties in the available information. Because our examples focus on the risk associated with the smaller asteroid size regime, the analysis is limited
to diameters 300 m and smaller. An upper diameter of 300 m was
selected to safely bound the impact risk at the smaller sizes and,
based on current research, avoid the uncertainties associated with
regional and global effects. The exclusion of regional and global effects is discussed further in Section 3.
2. Probabilistic asteroid impact risk model
The PAIR model presented is an extension of the work first
described in Motiwala et al. (2015). The current version of the
model has been extended to include an improved breakup model
that estimates burst altitudes from specific energy deposition rates
for a variety of fragmentation characteristics, a thermal radiation
damage model, and a more flexible scenario generation approach
that allows input parameters to be defined using a variety of
uncertainty distributions. Currently, the model execution consists
of three phases: impact scenario generation, impact consequence
modeling, and risk analysis of the combined hazard results. The
model generates sets of potential impact scenarios, performs physical modeling of the asteroid’s entry, breakup, and ground damage
for each case, computes the population affected by each modeled
impact, and then quantifies the aggregate risks by weighting the
relative damage probabilities with corresponding impact frequency
estimates.
2.1. Impact scenario generation
The scenario generator consists of a Monte Carlo framework
where distributions of input variables are defined and sampled to
produce a specified number of impact cases for analysis. The input variables include asteroid properties (such as size, density, and
strength), impact parameters (such as velocity, entry angle, and
impact location), and some modeling parameters (such as ablation
coefficient and luminous efficiency) that can be varied to capture
either modeling uncertainty or ranges of potential physical variations. Each parameter can be held constant, treated as uniformly
or normally distributed, or given a custom distribution.
The asteroid size can either be assigned directly through a
diameter distribution, or can be determined from sampled Hmagnitude and albedo values. As will be discussed in the breakup
modeling section, the asteroid strength is represented using an initial aerodynamic breakup pressure, along with a fragment strength
scaling parameter that controls the successive fragmentation rate
once breakup begins. The ablation coefficient varies the mass ablation rate used to compute atmospheric energy deposition in the
entry-modeling phase, and the luminous efficiency is used to estimate the potential thermal radiation.
Upon execution, the scenario generator samples the appropriate distributions and generates a list of the user-specified number
of impacts. The specific input distributions used for current assessment results are presented in Section 3 below.
2.2. Atmospheric entry and breakup modeling
For each impact scenario, the asteroid’s entry and breakup are
modeled to estimate an airburst altitude or ground impact based
on the energy deposited in the atmosphere through drag and ablation. The energy deposition is modeled using the fragment-cloud
model (FCM) from Wheeler et al. (2017), which combines a discrete, progressive fragmentation approach with releases of dispersing debris clouds. This approach is able to reproduce the types of
108 D.L. Mathias et al. / Icarus 289 (2017) 106–119
large flares observed in the Chelyabinsk meteor event and also allows for larger, stronger fragments to penetrate further and burst
lower. A summary of the modeling approach is included here, and
the Wheeler et al. (2017) and Register et al. (2017) references contain a complete discussion of the fragmentation and energy deposition models.
The entry flight dynamics are modeled by integrating the
single-body equations (Opik, 1958):
dm/dt = −0.5ρairv3Aσ (1a)
dv/dt = −0.5ρairv2ACD/m − gsinθ (1b)
dθ /dt = (v/(RE + h) − g/v)cosθ (1c)
dh/dt = vsinθ (1d)
where m is the mass of the asteroid, v is the velocity relative to
the atmosphere, θ is the flight path angle, h is the altitude, t is
time from the initial interface, g is acceleration due to gravity, ρair
is the local atmospheric density, RE is the radius of the earth, A
is the cross-sectional area of the object, CD is the drag coefficient,
and σ is the ablation coefficient.
The initial flight integration begins 100 km above the surface of
the earth and continues until the stagnation pressure exceeds the
object’s aerodynamic strength:
ρairv2 > aerodynamic strength (2)
Once the asteroid’s aerodynamic strength is exceeded, the parent body is broken into a given number of individual child fragments and an aggregate debris cloud. The model allows the number of fragments and cloud mass per split to be specified at run
time, with a baseline setting of two fragments and 50% cloud
mass for each fragmentation event. The child fragments are treated
themselves as intact bodies with an increased strength given by
Schild = Sparent(mparent/mchild )
α (3)
where S represents the strengths, m the masses, and α is the
power law strength-scaling exponent. The child fragments then
continue to fly until they reach their respective failure strengths
and again break into the prescribed number of fragments and
cloud mass. The process continues until the fragments reach the
ground or become slow or strong enough to no longer reach
the breakup condition. Following the approach of Hills and Goda
(1993), the debris cloud mass is treated aerodynamically as a single deformable body that spreads and slows under a common bow
shock. The lateral spread rate of the body’s effective drag area is
given by the expression (Hills and Goda, 1993):
vdispersion = vcloud (3.5ρair/ρcloud )
1/2 (4)
where vcloud is the velocity of the debris cloud, decelerating with
each integration step, and ρcloud is the density of the initial asteroid material. Passey and Melosh (1980) derive a similar expression for fragment spread rates, which differs only slightly in the
value of the constant being 3.0 (vs. 3.5) for analogous input assumptions. Laurence (2006) has also studied the relative motion of
fragments in hypersonic flight. The mechanics of the cloud debris
dispersion are important to the energy deposition and subsequent
ground damage that drives the risk. Sensitivity studies in Wheeler
et al. (2017) indicate that dispersion rate coefficients within the
cited 3–3.5 range or higher have little effect on energy deposition
(for the baseline assumption of 50% cloud mass), although reducing the coefficient by more than half could lower burst altitude estimates. A thorough discussion of the energy deposition modeling
can be found in Wheeler et al. (2017).
As the fragments and debris clouds decelerate and ablate according to the above equations, the total loss of kinetic energy of
Fig. 1. Energy deposition curve computed using the fragment-cloud model (solid
line), showing burst altitude at peak (dashed line). The case shown is for a 100 m
diameter asteroid with a density of 2.5 g/cc and aerodynamic breakup strength of
1 MPa entering at a velocity of 20 km/s and angle of 45° from horizontal.
all components is tracked and summed to estimate the energy deposited in the atmosphere per unit altitude. The point of maximum
energy deposition is taken to represent the approximate burst altitude of the impact case. Fig. 1 shows a sample FCM energy deposition curve and the burst altitude associated with the peak. If the
asteroid fails to break up or the energy deposition rate is still increasing at ground level, then the case is considered to be a ground
impact with a burst altitude of zero. The burst altitude is then used
to estimate the damage areas due to blast overpressure and thermal radiation as described below.
2.3. Blast overpressure damage
Blast waves from an asteroid airburst or ground impact can
cause varying levels of infrastructure damage or casualties over
large areas. The PAIR model uses curve fits of nuclear test data
(Hills and Goda, 1993) to estimate the blast damage radius for a
given burst energy and altitude. The burst altitude is based on the
peak energy deposition point as described in Section 2.2, and the
burst energy is assumed to be equal to the full initial kinetic energy of the asteroid. While in reality the energy going into the
blast will be some fraction of the deposited energy, it is not well
known what energy fraction is appropriate compared to the traditionally utilized static nuclear sources. Assessments presented in
Hills and Goda (1993) and Stokes et al. (2003) assumed that all
of the kinetic energy contributes to the blast for their calculations,
while the assessment in Toon et al. (1997) assumed a 50% contribution with the rest going into other energy modes such as radiation. Hills and Goda (1998) later suggested that an amount less
than 100% would be appropriate. However, Boslough and Crawford
(2008) and Aftosmis et al. (2016) have shown that a moving energy
source does not act identically to the static airburst assumed in the
nuclear scaling relations, and tends to produce a ground footprint
that appears stronger than an equivalent static source. For this reason, putting 100% of the kinetic energy into the blast scaling relations seems to correlate better with the high-fidelity simulations
and is used herein to provide a bounding, worst-case assessment.
Following previous blast damage conventions (Stokes et al.,
2003; Hills and Goda, 1993; Glasstone and Dolan, 1977), the baseline ground damage area is taken as the region within which
overpressure levels exceed 4 psi. The 4-psi ground damage radius,
Rground, is estimated using the scaling relation that Hills and Goda
(1993) derived from nuclear sources in Glasstone and Dolan (1977),
given by:
Rground = 2.09h − 0.449h2E−1/3 + 5.08E1/3 (5)
D.L. Mathias et al. / Icarus 289 (2017) 106–119 109
where E is the impact energy in megatons (Mt), and h is the burst
altitude in kilometers.
2.4. Thermal radiation damage
In addition to blast overpressures, large bursts or impacts can
release damaging levels of thermal radiation. The PAIR model
adapts the model from Collins et al. (2005), which estimates the
radius within which thermal radiation exceeds damage-causing
limits. The model computes these distances based on impact energy and a luminous efficiency parameter that defines what fraction of the energy is emitted as thermal radiation. Collins notes
that this parameter is poorly constrained, and cites a range of 1e−4
to 1e−2, with 0.003 taken as their assumed baseline value. For the
current assessment, the luminous efficiency is sampled as an uncertain parameter for each impact case, as described in Section 3.1.
The original formulation from Collins assumes that the thermal radiation emanates from a fireball plume of a ground impact, and
distributes the available radiation energy over the hemispherical
area above the ground to determine the exposure threshold radius:
r =

ηE/2πi , i = ZiE1/6 (6)
where r is the threshold radius from the burst or impact point,
η is luminous efficiency, E is the impact energy, and i is
the thermal exposure (total heating per unit area) associated
with a given damage severity. The damage severity is computed for a given impact energy based on energy scaling from
the corresponding thermal exposure threshold for a 1-megaton
explosion, Zi. The thermal exposure level corresponding to 3rd
degree burns (Zi = 0.42 MJ/m2/Mt1/6) is taken as the baseline
damage area threshold, with an option to output 1st degree
(Zi = 0.13 MJ/m2/Mt1/6) and 2nd degree (Zi = 0.25 MJ/m2/Mt1/6)
burn threshold areas as well. Implicit in the formulation is the
time over which the energy is radiated. As Collins points out, larger
explosions radiate more energy than smaller explosions, but the
time over which the energy is released is also longer. The i in
the denominator of the radius expression is the result of the exposure time. All results are scaled from the 1-megaton results derived
by Glasstone and Dolan (1977).
The PAIR model extends the formulation to include airburst
cases by taking the Collins hemispherical radius, r, centered at the
burst altitude, h, and computing the intersecting ground radius
as Rground =

r2 − h2 accordingly to determine the thermal damage area. Although thermal energy would realistically be emitted
spherically around the airburst, the current model maintains the
hemispherical energy distribution assumed in the original formulations, meaning the spherical exposure threshold radius, r, may
be overestimated by a factor of √
2 for airburst cases. This pessimistic convention was initially adopted to bound the thermal radiation risks and to remain consistent with the published Collins
model results for ground impact cases. Practically, this bounding
assumption does not significantly skew the results presented here,
as the blast overpressure is the dominant damage-causing source
for the majority of impacts in the size range considered. A quantitative discussion of the radiation surface assumption, as well as a
comparison of blast versus thermal radiation damage, is included
in Section 4.1.
2.5. Local affected population estimates
Blast and thermal radiation damage radii are computed for each
specific impact case, and the larger of the two is taken as the local
damage area. The damage area is centered at the coordinates sampled for the given impact, and the local population within that area
is counted as the affected population in the damage assessment.
Fig. 2. Notional diagram of how affected population is computed from the population grid. The figure on the left shows a damage circle over a 6 × 6 grid-cell region.
The right-hand figure shows a close-up of one of the boundary grid cells, divided
into a user-specified number of sub-cells, with the shaded overlay showing the subcell areas contributing to the population fraction counted within that cell.
Although actual damage severities and casualty rates would fall off
gradually with distance, including both survivors inside the defined
damage region and some casualties occurring beyond it, the 4-psi
and 3rd degree burn thresholds are assumed to provide a level of
severity at which much of the population would be substantially
affected by the damage or, alternately, provide a reasonable balance between survivors within and casualties without. Use of the
population within a distinct damage area also allows for additional
risk sensitivity analyses to be performed in the post-process phase
using localized population densities.
The local populations are computed based on gridded population data from the Socioeconomic Data and Applications Center
(CIESEN, 2016). The data set used in the current assessment provides population counts within 2.5-arc-minute grid cells, based on
UN-adjusted census data from the year 2000. Damage areas typically overlap portions of multiple grid cells. The fraction of the
population affected in each cell is computed by dividing the cell
into a user-specified number of sub-cells, and scaling the cell’s
population count by the fractional area of sub-cells with center
points inside the damage radius. Fig. 2 shows a diagram of a damage area over a grid region and an example of the sub-cell areas
used to compute the population fraction for a cell on the damage
boundary. This gridded population count is repeated for each of
the impact scenarios, using the specific computed damage area and
sampled impact coordinates for each case. Results for each simulated strike are then stored and processed in the final risk assessment phase.
2.6. Risk analysis
Impact consequences are converted into the related risk by
weighting the affected population results by the likelihood of such
impacts. Because the current assessment focuses on the ensemble
risk rather than a specific impact scenario, the frequency of impact
is determined using the general asteroid population estimates published by Harris (2015). The results provide bounding estimates on
the risk, including the uncertainty around the lack of complete asteroid population knowledge, and can be used to identify how the
risk uncertainty changes based on additional information about the
population.
3. Ensemble asteroid threat assessment of sub-300 m impactors
We applied the PAIR analysis approach to a statistical range
of hypothetical impactors up to 300 m in diameter to establish
bounds on the resulting risk. For this assessment, risk is defined
via the affected population within the damage radii as described
above.
110 D.L. Mathias et al. / Icarus 289 (2017) 106–119
The current discussion focuses on demonstrating the PAIR
model framework and illustrating how the model can be used to
provide insight into specific decisions—in this case, considering the
minimum size of objects that constitute a threat. Since the focus
is on the small object size range, we restrict the assessment to
objects 300 m in diameter or smaller and only consider the local
damage effects associated with blast overpressure and thermal radiation. Reinhardt et al. (2016) compare blast, fireball, thermal radiation, crater size, and seismic impact severity for a range of impactor sizes. Blast and thermal radiation bound the damage radius
across the size range considered here. Stokes et al. (2003) include
“global effects” to represent widespread consequences, for example from impact ejecta, but these are limited to impactors above
1 km in size. Stokes et al. (2003) also include impact tsunami consequences, based on Chesley and Ward (2006), which are maximum for impactors of 300 m. However, results from the 2016
NASA-NOAA Asteroid-Generated Tsunami (AGT) Workshop1 indicate that ensemble impact tsunami risk is lower than the corresponding blast/thermal risk for objects 300 m and below (Morrison
and Venkatapathy, 2017). The local airburst and thermal radiation
results shown in this section bound the total risk for the size range
considered. Inclusion of the tsunami and global effects models is
left as future work.
We begin by discussing the uncertainty distributions used for
each sampled input parameter, which were all defined based on
established ranges from the literature. Model convergence is then
covered to establish the number of Monte Carlo samples needed
to ensure that the uncertainty in the output is driven by the input
distributions and not a bias from insufficient sampling. Results are
shown using the average expected casualty metrics and compared
with the full risk probability distributions enabled by the PAIR assessment approach. We then consider the question of how small of
an asteroid constitutes a significant threat to humans. This requires
the introduction of a specific risk tolerance posture, and we use
the model results to answer the question in that context. Lastly,
we revisit the size-of-threat question considering the uncertainty
when the object size is based only on optical wavelength photometry. The results illustrate how the current risk model can provide
information highlighting the sensitivity of the results to the underlying assumptions and potentially identify areas where increased
knowledge could improve our understanding of the threat.
3.1. Input parameter distributions
For this assessment, 30 million stochastic impacts were modeled. Each available input parameter was allowed to vary freely
from case to case, sampled randomly from a baseline set of uncertainty distributions derived from the literature. The distributions
and ranges are described below and shown in Fig. 3.
Size: Asteroid sizes were determined based on sampled Hmagnitude and albedo values. The H-magnitude was sampled uniformly from 20 to 30, and the albedo was sampled from the
NEOWISE distribution (Mainzer et al., 2011). The diameter was
then computed from the sampled values as: D = (1.326 × 106) ×
10−H/5/
√pv, where D is in meters, H is the absolute magnitude,
and pv is the albedo (Bowell et al., 1988). Sampling in H-magnitude
space naturally skews the size distribution toward the smaller end
of the range. This helps the convergence of the results as will be
discussed in Section 3.2 and is appropriate for the specific applications discussed in Sections 3.4 and 3.5.
Density: Current estimates of asteroid densities remain notably
uncertain and debatable, since measurements of pre-entry asteroid
1 Second International Workshop on Asteroid Threat Assessment: Asteroidgenerated Tsunami (AGT) and Associated Risk Assessment, Seattle, WA, August 23-
24, 2016. Available: https://tsunami-workshop.arc.nasa.gov/workshop2016/.
Table 1
Base meteorite densities (Macke, 2010) and fall fraction (Spohn et al., 2014)
used to derive asteroid density distribution.
Meteorite Type Meteorite Density (g/cm3) Fraction (%)
H 3.38 ± 0.19 34.9
L 3.30 ± 0.12 38.9
LL 3.19 ± 0.14 9.3
CM 2.27 ± 0.13 4.3
IAB 6.75 ± 1.84 2.5
IIIAB 7.15 ± 0.57 2.4
EUC 2.84 ± 0.13 3.4
DIO 3.12 ± 0.19 1.1
HOW 2.86 ± 0.11 3.4
densities are extremely difficult and meteorites that survive to the
ground are not necessarily representative of asteroids. The distribution of available measured asteroid densities is asymmetrical and
sensitive to selection effects. However, the distribution of macroporosities derived from asteroid density measurements is symmetrical and relatively insensitive to selection choices. Our nominal
input density distribution is based on the density distribution of
meteorites, modified by a likely macroporosity distribution. The
meteorite density distribution was calculated based on measured
density distributions of different meteorite classes (Macke, 2010)
and the relative fraction of each class by falls (Spohn et al., 2014),
as shown in Table 1. We generated a base meteorite density distribution by defining a Gaussian distribution for each meteorite
type based on its mean density and standard deviation, and subsequently combining the Gaussians in the ratios specified by the
fall fractions. This base meteorite density distribution was modified by an asteroid macroporosity distribution based on the values
tabulated by Carry (2012). Based on the 97 asteroids having density measurements with relative uncertainties ≤ ± 50% and macroporosities greater than zero, we fit a Gaussian distribution with a
mean of 34% and sigma of 18% to model the macroporosity. Additional limits were placed on the final asteroid density distribution requiring that the asteroid densities must fall between 1.1 and
7.5 g/cm3 and that the macroporosity must be within the range of
1–70%.
Strength: The initial aerodynamic breakup strength was sampled from a logarithmic distribution from 0.1 to 10 MPa (Popova
et al., 2011). Published estimates of fragment strength scaling exponent values vary widely from as low as 0.03 (Yoshinaka et al.,
2008) to 0.57 (Cotto-Figueroa et al., 2016). Based on the summary
by Popova et al. (2011), a uniform distribution between 0.1 and 0.3
was used for the current assessment. Sensitivity studies of FCM
fragmentation parameters have also shown these ranges to produce reasonable variations of energy deposition curves to represent uncertainties in breakup mechanics due to asteroid composition (Wheeler et al., 2017).
Impact velocity: This assessment uses an impact velocity distribution derived by Greenstreet et al. (2012). Greenstreet et al. simulated near-Earth asteroid populations by modeling the injection of
objects into the near-Earth region due to resonances, utilizing an
approach very similar to Bottke et al. (2002). Due to increases in
computing power, the Greenstreet model was able to simulate 6
times more objects with finer time-step resolution, permitting improved statistics and improved modeling of high-speed planetary
encounters.
Entry angle: Entry angles are sampled from a cosine distribution that weights entry angle probabilities toward 45°. Studies have
shown 45° to be the most probable entry angle, with vertical or
very shallow entries much less likely (Shoemaker, 1961). The distribution is given by θ ◦ = (90◦/π)cos−1(2U − 1), where U is a random number sampled uniformly between 0 and 1.
Impact location: For ensemble impact assessments, the model
randomly selects longitude and latitude so that the entire area of
D.L. Mathias et al. / Icarus 289 (2017) 106–119 111
Fig. 3. Histograms of the input parameter distributions sampled to produce the impact cases for the presented risk assessment.
the globe is evenly sampled. Gallant et al. (2009) used the debiased NEO population model of Bottke et al. (2002) to examine the
spatial distribution of objects impacting earth and found that, averaged over the year, there was less than a 5% latitudinal variation
after correcting for even areal coverage.
Ablation parameter: The ablation coefficient, σ, used in the
flight integration was sampled from a logarithmic distribution between 3.5e−10 and 7e−8 kg/J. This broad range was selected to
bound the uncertainty in this parameter by combining the values from the asteroid impact literature (Hills and Goda, 1993,
1998; Collins et al., 2005; Shoemaker, 1961) with assessment of
the material heat of ablation (Chen, 2016) and the uncertainty
in the aerothermodynamic heating environment (Johnston et al.,
1995).
Luminous efficiency: The luminous efficiency parameter used in
the thermal radiation model was also sampled logarithmically from
a large, highly uncertain range found in the literature. A range from
3e−4 to 3e−2 was selected based on the 1e−4 to 1e−2 range from
Collins et al. (2005), but shifted slightly to center logarithmically
on the Collins et al. nominal 0.003 value.
3.2. Convergence
Given the input distributions, sufficient scenarios must be run
so that the resulting statistical distributions are meaningful, covering the full range of potential outcomes thoroughly enough that
the results do not continue to change with additional cases. For
the current assessment, 30 million impact cases were run, with
112 D.L. Mathias et al. / Icarus 289 (2017) 106–119
Fig. 4. Convergence of the full impact set with all asteroid sizes. The plot shows the
percentage difference between a running average of the location-specific gridded
population counts and the average population within the average damage area, as
a function of number of impact cases.
varied asteroid properties, entry parameters, impact locations, and
model parameters all sampled from the distributions described in
the previous section. The various damage areas produced by the
range of asteroid parameters sampled converges quickly to an average value within a fairly small number of cases. The location of
impact, however, causes the results to vary significantly; impacts in
the ocean infrequently cause damage while similar impacts in large
metropolitan regions result in large affected populations. Therefore, the affected population count serves as the best convergence
measure at hand. If too few points are sampled, then a single strike
over a densely populated area could shift the aggregate results, indicating that the solution has not yet converged. To establish convergence for a given size range, the local affected population results produced by the model are compared to results assuming a
uniformly distributed population. Specifically, a running average of
the location-specific gridded population counts is compared with
a running average of the damage areas multiplied by the average
world population density, as a function of the number of impact
cases.
Fig. 4 shows the population convergence trends for 30 million
simulated impacts. The overall results are within ∼5% of the average values by 1 million realizations, and are within a fraction of a
percent by the full 30 million realizations. Fig. 5 shows the average
population convergence for four representative size bins within the
set: 40–60 m, 80–100 m, 180–200 m, and 280–300 m.
As expected, the smallest objects require the most scenarios to
converge because the resulting damage areas are the smallest, differ the most based on impact location, and require more points to
cover enough overall surface area around the globe. The 40–60 m
bin, which represents the smallest sizes able to produce statistically notable damage, is within 6% of the average values by the 2.7
million realizations modeled in that size range, while the larger
bins are all within ∼1% by a few hundred thousand realizations.
Overall, these results are sufficiently converged that further variations from additional cases do not noticeably change the ensemble
risk results presented.
3.3. Local damage source results
Before considering the ensemble risk results, we investigate the
relative contributions of the local damage sources and the sensitivity of the expected casualty estimates to thermal radiation modeling assumptions.
Dominant damage source: In the beginning of Section 3 we
posited that, for the size range presented, the impact consequences
are bounded by blast and thermal radiation effects. We follow
Stokes et al. (2003), Hills and Goda (1993), and Reinhardt et al.
(2016) in the use of 4-psi overpressure as the criterion for establishing casualties. Boslough and Crawford (2008) use blast-induced
winds as a metric for damage assessment but, as outlined in
Glasstone and Dolan (1977), the wind speed and overpressure both
relate to the source energy. Fig. 6 shows the fractional contribution
of blast and thermal radiation hazards in terms of the percentage
of the Monte Carlo scenarios as a function of size. In other words,
the plot shows the fraction of cases, within a given size range,
for which blast or thermal is the dominant damage source (or for
which neither hazard produce any damage). Below 50 m, most of
the simulated scenarios cause no casualties, and above 50 m blast
damage is the dominant source. As noted by Collins et al. (2005),
the luminous efficiency parameter in the thermal radiation model
is highly uncertain, and the current results include the large variation discussed in Section 3.1. If the nominal value of 0.003 is used,
blast damage exceeds thermal radiation damage for every scenario
across the sizes presented.
Thermal radiation: Thermal radiation damage is defined by
the level that would cause 3rd degree burns, as described in
Section 2.4. In that discussion, an assumption that the thermal energy is distributed over a hemispherical surface was made. Strictly
speaking, this is only true for ground impacts, and bursts high
above the ground would distribute their energy spherically. However, it is not clear a priori if a specific case is better represented
by a spherical or hemispherical distribution area—in reality most
fall somewhere in between. We use the hemisphere for all cases
because it serves as a “worst case.” Fig. 7 shows a comparison of
the expected casualty results using the spherical vs. baseline hemispherical assumption. The factor of two in radiated energy flux
makes less than 10% difference in the cumulative results. This is
because the thermal radiation is overshadowed by the blast damage for the current size range.
3.4. Ensemble risk assessment results
A traditional approach to quantifying impact risk involves computing an average expected casualty rate (Ec) for impacts from asteroids up to a given size. In this case, the casualty results from
impacts within a given size range are simply averaged and then
multiplied by the expected impact frequency for objects of that
size. Fig. 8 shows an example of the cumulative and differential
expected casualty curves generated from the current impact risk
assessment, taking the affected population estimates as the “casualties”. The casualty rates plotted on the left represent the cumulative risk from all impacts up to the given size, while the right plot
represents the differential casualties per meter of impactor size.
While such an approach seems to provide a very simple cutand-dry metric for evaluating aggregate risk, it fails to represent the true breadth of risk potential for low-probability, highconsequence events like asteroid strikes. Figs. 9 and 10 show the
variation in possible outcomes for four asteroid size ranges, assuming an impact occurrence within each size group. These conditional
strike results reveal the variation in consequences produced by including uncertain distributions of asteroid impact parameters, prior
to accounting for the overall probability of such asteroids striking Earth. Fig. 9 shows the relative likelihoods of different affected
population ranges for each size. Even when a strike is assumed to
occur, the results exhibit a bimodal behavior. The most likely result for any of these sizes is that no population will be affected.
However, should the strike occur with consequences, those conse-
D.L. Mathias et al. / Icarus 289 (2017) 106–119 113
Fig. 5. Convergence of local affected population for four size bins within the 30 M impact cases. Plots show the percentage difference between a running average of the
location-specific gridded population counts and the average population within the average damage area, as a function of number of impact cases. The number of impact
cases plotted reflects only the cases sampled within the given size range rather than the total cases in the overall set. The dashed lines show the final value that each bin
has converged to over all of cases sampled in this assessment.
quences can become quite severe. For a given impactor size, the
severity of the hazard depends significantly on the impact location. Fig. 10 compares the min, mean, and max affected population results for the same size ranges, using the local population
at specific impact sites vs. assuming an average world population
density within the damage areas. While the means naturally remain the same, the maximum consequence levels are an order
of magnitude higher when the potential to hit highly populated
areas is included. Similarly, even 300-m objects that always produce relatively large damage areas maintain a high probability of
striking the open ocean and may cause no casualties. These important distinctions are lost when distilling results to mean values or using single, representative strike parameters. By using the
mean, the range of possible outcomes is completely masked and
the full range of information is not represented in an actionable
manner.
While Figs. 9 and 10 move beyond the simple average outcome,
a more complete view of the total risk picture is desired. To this
end, Fig. 11 shows a contour plot of the likelihood of exceeding
various damage levels assuming an impact occurs within each size
bin considered. As in Fig. 9, these conditional probabilities reflect
the range of outcomes stemming from the input parameter distributions and do not account for the frequency with which an
impact in that size range may occur. However, rather than giving
the relative probabilities of each individual damage range, the exceedance probabilities give the complementary cumulative probability of affecting at least the given population threshold or greater.
The use of cumulative probabilities minimizes the arbitrary dependence of the probability values on the width or number of bins selected and instead enables clearer evaluation of meaningful thresholds. Instead of giving a single average value, Fig. 11 shows the
full range of potential consequences, and their relative likelihoods,
given the uncertain variations of asteroid properties and entry parameters.
To assess the absolute level of risk posed by different asteroid
sizes, the conditional strike results are then weighted by the likelihood of such impacts occurring. To do this, the estimated annual
impact frequency for a given size range, interpolated for each size
bin from the values published in Harris (2015), is multiplied by the
corresponding horizontal slice through Fig. 11. In order to avoid the
bin-width-dependence that comes with assigning impact probabilities for a given size range, the probabilities are now also presented
as cumulative over size as well as damage, representing the probability per year of an impact of a given size or smaller affecting a
given population number or higher. The resulting cumulative annual damage exceedance contours are shown in Fig. 12.
114 D.L. Mathias et al. / Icarus 289 (2017) 106–119
Fig. 6. Dominant damage source as a fraction of Monte Carlo scenarios for asteroid
size of up to 300 m.
Using these probability contours, one can then add a risk posture in order to address the question of what size constitutes a
considerable threat. In practice, it can be challenging to define a
distinct, quantitative risk posture, but for this example it is assumed that the threshold for constituting a threat is a one-ina-million annual probability of affecting 10,000 people or more.
The dashed arrows in Fig. 12 illustrate how the results can provide a size threshold for such a risk tolerance. The bold black
line shows the one-in-a-million annual probability contour, and
the intersecting dashed lines show the size at which the sample
10,000-person threshold exceeds that probability. Based on these
assumptions, objects of ∼65 m or smaller would not exceed the
given risk threshold, while larger objects would be considered a
threat.
In these ways, the cumulative annual damage exceedance figures produced by the PAIR model can provide a full representation
of knowledge, or the limits of knowledge, in terms of the uncertainty distribution, while also enabling quantitative answers to be
distilled for specific questions.
3.5. Effects of observational size uncertainty on threat evaluation
This approach to representing risk uncertainties can further
be used to evaluate how different parameter distributions, correFig. 7. Comparison of cumulative expected casualties as a function of asteroid size
for spherical and hemispherical thermal radiation surfaces. The results include both
blast and thermal radiation damage.
sponding to alternate states of knowledge, can influence the answer to the hypothetical size threshold question. Such results can
be used to prioritize investments to gain additional knowledge in
the context of a threat, or can help characterize the likely limits
to the benefit of gaining additional data. As an example, we consider a case study looking at how the results change given imperfect knowledge of object diameter.
To investigate the effects of uncertain size estimates based
solely on optical observations, Fig. 13 shows the same simulation results from Fig. 12 re-binned by the optically observable Hmagnitude rather than by diameter. These results inherently incorporate the size uncertainty associated with each H-magnitude bin,
because the scenario generation process produces ranges of actual asteroid diameters based on the sampled albedo distribution.
When the same sample risk posture is applied to H-magnitudebased results in Fig. 13, an H-magnitude threshold of around 26.5
is obtained. Converting this to an equivalent size threshold using
the conventional average albedo value of 0.14 yields a threshold of
less than 20 m, rather than the 65 m found from the actual impact
diameters shown in Fig. 12. Conversely, using the average albedo
to convert the 65 m size threshold to optical survey criteria yields
an H-magnitude threshold of ∼24, rather than the ∼26.5 threshold obtained from Fig. 13. Put another way, the results show that
the uncertainty involved with inferring an object diameter based
Fig. 8. Expected casualty results (average casualties per year). On the left, the casualties are presented as cumulative values for impactors up to the given diameter threshold.
On the right are the differential expected casualties per unit size (casualties per year per meter). The right plot is the slope of the curve on the left.
D.L. Mathias et al. / Icarus 289 (2017) 106–119 115
Fig. 9. Conditional probabilities of impacts within four size ranges causing varying
levels of population damage, given an impact in each size group. The bars represent
the relative probabilities of each population range, while the lines represent the
cumulative probability up to each range.
Fig. 10. Comparison of the minimum, mean, and maximum affected population
consequences given impacts within each size range, assuming either an average
world population density (dashed lines) or using local populations at specific impact sites (solid lines).
Fig. 11. Probability (color contours) that an impact from an asteroid within a given
10-m size range (y-axis) will affect at least a given number of people or more (xaxis), assuming that an impact of that size occurs.
Fig. 12. Probability per year (color contours) of an asteroid up to a given size
threshold (y-axis) impacting Earth and causing damage that affects at least a given
population threshold or more (x-axis). The black lines show the 10−5, 10−6, and
10−7 probability contours. The dashed arrows show an example of how impact risk
probability distributions can provide a size threshold for a sample risk tolerance
level of a one-in-a-million probability per year of an impact affecting at least 10,000
people or more. For the modeling and parameter assumptions used in this assessment, the size threshold for this sample risk posture is around 65 m in diameter.
Fig. 13. Probability per year (color contours) of an asteroid up to a given absolute
brightness (y-axis) impacting Earth and causing damage that affects at least a given
population threshold or more (x-axis). Absolute magnitudes were converted to object diameters using the full distribution of albedos discussed in § 3.1. The bold
black line shows the one-in-a-million probability contour used in the sample risk
posture presented in Fig. 12.
solely on an optical detection lowers the definition of the sample
threat threshold from 65 m to an equivalent size threshold below
20 m. An optical system of sizing objects would have to search for
objects down to a 20-m equivalent size to capture the same risk
exposure as a “perfect” sizing system would by detecting objects
of 65 m.
This example demonstrates how the PAIR approach of using full
probability distributions enables key uncertainties to be appropriately incorporated into more meaningful and informative risk metrics. Use of mean parameter values and mean risk metrics, on the
other hand, can neglect or misrepresent the effects of such uncertainties on potential mitigation planning or survey decisions.
4. Discussion
In this section we discuss the current results in the context of
specific, previous efforts. First, the ensemble risk results are compared to those from Stokes et al. (2003). This comparison illus-
116 D.L. Mathias et al. / Icarus 289 (2017) 106–119
Fig. 14. Comparison of cumulative expected casualty results from the current PAIR assessment and the Stokes et al. (2003) assessment, showing the cumulative contributions
of key modeling and input differences incorporated into the PAIR assessment. Each curve, from top to bottom, shows the results with the addition of the variation listed in
the legend, while also maintaining all the previous variations added in the lines above. The solid black line near the bottom is the final PAIR assessment result, including all
modeling updates and input parameter distributions.
trates the differences through incremental changes in modeling assumptions and input values. Additionally, this build-up will highlight the modeling aspects that most significantly drive the results.
We then discuss the airburst altitude predictions as applicable to
the risk assessment. The Tunguska event is used as a sample case
for comparing the current results to previous work. Throughout the
discussion we highlight opportunities for advancement in modeling, and conclude with a summary of our planned future work.
4.1. Comparison of risk results with previous studies
Of the existing quantitative impact risk assessment studies,
Reinhardt et al. (2016) also consider a stochastic asteroid impact population, but the results are presented for a different
size range of interest and are not directly comparable. Rumpf
et al. (2016) consider the specific potentially hazardous asteroid
(PHA) population and is also not directly comparable. Stokes et al.
(2003) provide the best opportunity for comparison to the current
results. Fig. 14 shows the average expected casualty results from
Stokes et al. (2003) as a function of cumulative asteroid size, compared with PAIR model results that incrementally incorporate the
key differences between the two assessments: impact frequency
updates (applied to distributed sizes), fragmentation modeling, and
distributions of asteroid and impact parameters rather than representative averages.
For validation, we first ran the PAIR model using the same input values and pancake breakup model assumptions used in the
Stokes assessment and reproduced the results, as shown by the
top dotted curve in Fig. 14. The subsequent curves represent incremental changes in the results from the previous curve (from
top to bottom). The second curve on the plot illustrates the change
in PHA population and corresponding impact frequency estimates
since 2003—this constitutes the single biggest increment between
Stokes et al. (2003) and the current work. In this result set, the
impact frequencies from Harris (2015) were used in conjunction
with distributed asteroid sizes as is done in the current PAIR approach. Rather than using the fixed size bins and frequencies, the
distributed sizes are binned by a chosen increment, and the impact frequency for each size range is interpolated from the Harris
values.
The next curve continues to use the Harris (2015) impact frequencies and adds in the FCM breakup modeling approach, with
50% cloud mass released in each break and a fixed fragment
strength scaling parameter of 0.2. Adding the discrete fragments
to the breakup models does increase the breakup altitude in some
cases, as discussed in Section 4.2 below, and shifts the resulting
casualties lower at the lower end of the size range. This difference becomes small in the cumulative results by 300 m. Compared
to other modeling changes, the difference between FCM and pure
pancake models is relatively minor.
In the next line, the aerodynamic strength and fragment
strength scaling distributions are introduced, compared to the
fixed nominal values in the previous curves. This addition reduces
the average casualties by nearly 50% across the range, including a
corresponding shift in the minimum damage-causing size. The next
curve introduces the density distributions described in Section 3.1.
This is a relatively small shift in results compared to the previous curve. The final curve includes the addition of all the other
input parameter distributions used in the presented PAIR assessment (i.e., velocity, entry angle, ablation coefficient, and luminous
efficiency).
The relative contribution of the individual parameters presented
does depend on the order in which the changes are applied. However, the build up of the differences in results between the current assessment and that of Stokes et al. (2003) is represented:
the largest difference is due to changes in impact frequency estimates; the changes due to the FCM modeling approach are relatively small; the biggest physical parameter influence is due to
density and strength; and all of the other distributions combine
for a relatively small increment for the current study. As discussed
in the previous section, the current analysis uses input parameters
based on the current literature. While these are not unique and alternate input formulations can shift the overall results, the trends
are likely to hold. A complete sensitivity study is currently underway and is part of the future work.
4.2. Energy deposition and airburst altitude
Energy deposition results depend both upon the breakup modeling assumptions employed and the inputs used to represent the
D.L. Mathias et al. / Icarus 289 (2017) 106–119 117
Table 2
Comparison of FCM airburst altitudes and ground damage radii for the Tunguska-scale (15 Mt) stone-type asteroid impact cases presented in
Chyba et al. (1993). The airburst altitude is determined from the modeled energy deposition peak, and the ground damage radius is the 4-psi
blast overpressure radius (computed from Eq. (5)). The all-cloud FCM results were computed using a cloud mass fraction input of 100%, and
the baseline FCM results were computed using 50% cloud mass and a fragment strength scaling α of 0.2 (representative of the values used in
the current PAIR assessment).
Entry Angle Airburst Altitude (km) Ground Damage Radius (km)
Chyba FCM all cloud FCM baseline Chyba FCM all cloud FCM baseline
15° 15 13.6 14.3 2.9 7.3 5.2
30° 11.5 10.8 12.6 12.5 13.9 10
45° 9 8.9 11.8 16.6 16.7 11.9
90° 6 6.5 10.5 18.5 18.4 14.4
asteroid characteristics. The resulting energy deposition curves can
be compared with observed lights curves for validation, and can be
used as energy sourcing into higher-fidelity blast modeling efforts.
However, for the ensemble risk model, the energy deposition curve
is used to determine the airburst altitude for the blast and thermal
radiation ground damage models.
Register et al. (2017) surveyed the primary energy deposition
modeling approaches based on the single-body entry equations
used in asteroid impact risk assessments. A finding from Register
et al. (2017) is that the “pancake” models can create representative flares but discrete fragmentation models are required to produce the multiple flares and varied energy deposition rates that are
common in observed meteor light curves. Wheeler et al. (2017) extended this finding to create the fragment-cloud model (FCM) used
in the current study. In Wheeler et al. (2017), the FCM was used
to reproduce the light curve for the 2013 Chelyabinsk meteor, and
achieved excellent matches through the incorporation of clouds
(“pancake” effects) and discrete fragments. This not only enables
the matching of multiple light curve flares and better matches the
peak deposition rate, but also moves beyond the assumption that
all asteroids are monolithic bodies prior to entry.
The inclusion of the fragments and clouds does change the
airburst altitude prediction when compared to the previous pancake models. To illustrate the effects of these differences in our
risk model, we compared FCM results with the assessment of
Tunguska-class stony asteroid airbursts performed by Chyba et al.
(1993). Table 2 shows the airburst altitude results for a 58 m diameter, 5.6e8 kg asteroid entering at 15 km/s (15 Mt of impact energy) for four entry angles (15, 30, 45, 90°). For the spherical bolide
assumption used in the FCM, a density of 5.48 g/cc was needed to
obtain an equivalent mass and energy as the 3.5 g/cc used with the
cylindrical bolide assumption in the Chyba et al. assessment. Results from the current FCM reasonably match those of the pancakestyle Chyba model when the breakup is modeled as a single cloud,
using the−higher breakup pressure (23.5 MPa) assumptions employed by Chyba et al.2 The airburst altitudes increase with the
addition of discrete fragments and multiple, smaller cloud masses,
which are required to match most observed light curves. These results do shift compared to the Tunguska literature, but that is not
surprising since the existing literature is largely based on similar
modeling assumptions, particularly the assumption that the object
was monolithic and disrupted in a single flare.
Table 2 also compares the resulting 4-psi ground damage radius corresponding to the different airburst altitudes. Again, when
using similar assumptions, the damage ranges agree quite well except at the very shallow entry angle. In this case, the relatively
small change in airburst altitude translates into around double the
ground damage area. The ground damage from a low 15 Mt burst
2 The “pancake” portion of the current model is derived from the work of Hills
and Goda (1993) and the results reproduce their results precisely when utilizing the
same inputs.
can be fairly sensitive to small changes in burst altitude, ranging
from maximum damage at the optimal burst height of 6 km, to
no damage by 16 km. It is also possible to see the opposite results where moderate changes in burst altitude cause very little
change in the ground damage. This phenomena results from the
non-linear shock propagation physics and does not specifically relate to the model choices in this study. Factors such as these highlight a sensitivity of the results to uncertainty in the knowledge of
asteroid compositions or structures. These sensitivities and uncertainties further support the use of result distributions over average values when assessing the range of consequences. The current
model permits a range of possible internal structures and makes it
possible to incorporate such aspects into the risk assessment. Work
is ongoing to compare the FCM energy deposition results to those
from detailed hydrocode simulations. The goal of this effort is to
relate physically meaningful strength and structure parameters to
FCM inputs so that a clear correlation can be made between asteroid characteristics and impact risk in a statistical sense.
4.3. Future work
The existing PAIR model framework is set up to modularly incorporate additional impact hazard models. Next steps include vetting the existing tsunami model against high-fidelity simulations.
Additionally, the global effects model of Stokes et al. (2003) will
be revisited to explicitly consider impact ejecta and the subsequent
effects. A formal study of the sensitivity of the risk results to the
modeling and asteroid parameter uncertainty is currently underway. This work will provide additional insight into the most important modeling parameters as well as providing a risk perspective
on the relative value of information in characterizing actual impact
scenarios. Finally, the relationship between the FCM aerodynamic
strength parameters and engineering measurements of strength is
being studied with hydrocode analyses.
5. Conclusions
A new probabilistic asteroid impact risk model was developed
to assess the risk that potential asteroid strikes pose to Earth’s
population. The model consists of a Monte Carlo framework that
simulates sets of impact scenarios sampled from distributions describing asteroid and entry characteristics. Such results can be used
to bound potential hazard levels for a wide variety of impact scenarios, assess overall ensemble risk from a general asteroid population, and provide insight into what properties or assumptions
most significantly affect the potential consequences or risk metrics
of interest.
In the current study, 30 million simulated impacts were used to
produce a nominal set of ensemble risk results, focusing on objects
of 300 m in diameter or smaller. The results were considered in
the context of identifying the minimum asteroid size that constitutes a threat for the given input assumptions. A hypothetical risk
118 D.L. Mathias et al. / Icarus 289 (2017) 106–119
tolerance level of a one-in-a-million likelihood per year of affecting at least 10,000 people was used to define a minimum impactor
size that corresponds to the threat level. For the baseline assumption set in the current analysis, an asteroid diameter of 65 m was
the smallest size to constitute that level of threat. When the ensemble of impact scenarios considered the effect of albedo uncertainty, optical search schemes would have to locate objects down
to H-magnitude of approximately 24, equivalent to a 20 m asteroid
assuming an average albedo of 0.14.
The primary impact hazards were shown for impactors up to
300 m. At this scale, seismic, tsunami, and global effects are not
important to the ensemble risk. Blast overpressure comprised the
dominant damage source and thermal radiation only contributed
for high values of luminous efficiency. Airburst altitude was also
considered using the Tunguska event as an example case. The current model reproduces the results from other researchers when
using the same inputs. However, the current model includes additional breakup modeling features that allow discrete fragments
and clouds to result from successive atmospheric breakup events.
This addition changes the burst altitudes slightly but also serves
as a proxy for a range of pre-entry asteroid structures (not requiring all pre-entry objects to be homogeneous monoliths). Such results reinforce the need to better understand asteroid composition
and structure in refining impact consequence assessments. Finally,
the results were compared to those from the 2003 NEO Science
Definition Team (Stokes et al., 2003), and a series of curves were
used to highlight the differences. The largest difference stems from
a change in the asteroid impact frequency estimates. Assigning distributions for asteroid density and strength, compared to fixed,
nominal values, made the next largest differences. Other modeling and characteristic combined for small relative differences from
the previous study.
Acknowledgements
This work was funded by NASA’s Planetary Defense Coordination Office (PDCO). The Gridded Population of the Word, Version 3 data were developed by the Center for International Earth
Science Information Network (CIESIN), Columbia University and
were obtained from the NASA Socioeconomic Data and Applications Center (SEDAC).
References
Aftosmis, M., et al., 2016. Numerical simulation of bolide entry with ground footprint prediction. 54th AIAA Aerospace Sciences Meeting.
Alvarez, L., Alvarez, W., Asaro, F., Michel, H., 1980. Extraterrestrial cause for the cretaceous-tertiary extinction. Science 208, 1095–1108.
Barringer, D., 1909. Meteor Crater (Formerly called Coon Mountain or Coon Butte) in
Northern Central Arizona, November 16. National Academy of Sciences, Princeton University Paper presented at the.
Boslough, M., Crawford, D., 2007. Low-altitude airbursts and the impact threat. In:
Hypervelocity Impact Proceedings of the 2007 Symposium HVIS.
Bottke, W., et al., 2002. Debiased orbital and absolute magnitude distribution of the
near-Earth objects. Icarus 156 (2), 399–433.
Bowell, E., 1988. Application of photometric models to asteroids. In: Asteroids II;
Proceedings of the Conference. Tuscon, pp. 524–556.
Carry, B., 2012. Density of asteroids. Planet. Space Sci. 73, 98–118 December. Provided by the SAO/NASA Astrophysics Data System.
Center for International Earth Science Information Network - CIESIN - Columbia
University. (2016) Socioeconomic Data and Applications Center (SEDAC).
[Online]. HYPERLINK "http://dx.doi.org/10.7927/H4639MP" http://dx.doi.org/10.
7927/H4639MP.
Chapman, C., Morrison, D., 1994. Impacts on the Earth by asteroids and comets:
assessing the hazard. Nature 367, 33–39.
Chapman, C., 2004. The hazard of near-Earth asteroid impacts on earth. Earth
Planet. Sci. Lett. 222 (1), 1–15.
Chen, Y., 2016. Thermal ablation modeling for silicate materials. 54th AIAA
Aerospace Sciences Meeting, AIAA SciTech Forum January[Online] Available:
http://dx.doi.org/10.2514/6.2016-1514.
Chesley, S., Chodas, P., Milani, A., Valsecchi, G., Yeomans, D., 2002. Quantifying the
risk posed by potential Earth impacts. Icarus 159 (2), 423–432.
Chesley, S., Ward, S., 2006. A quantitative assessment of the human and economic
hazard from impact-generated tsunami. Nat. Hazards 38, 355–374. doi:10.1007/
s11069-005-1921-y.
Chyba, C., Thomas, P., Zahnle, K., 1993. The 1908 Tunguska explosion: atmospheric
disruption of a stony asteroid. Nature 361, 40–44.
Collins, G., Melosh, H., Marcus, R., 2005. Earth impact effects program: a web-based
computer program for calculating the regional environmental consequences of
a meteoroid impact on Earth. Meteorit. Planet. Sci. 40 (6), 817–840.
Cotto-Figueroa, D., et al., 2016. Scale-dependent measurements of meteorite
strength: implications for asteroid fragmentation. Icarus 277, 73–77.
Ericson, C.II, 2005. Hazard Analysis Techniques for System Safety. John Wiley & Sons,
Inc., Fredericksburg, Virginia.
Gallant, J., Gladman, B., Cuk, M., 2009. Current bombardment of the Earth-Moon
system: emphasis on cratering asymmetries. Icarus 202 (2), 371–382.
Gersonde, R., et al., 1997. Geological record and reconstruction of the late
Pliocene impact of the Eltanin asteroid in the southern ocean. Nature 390,
357–363.
Glasstone, S., Dolan, P., 1977. The Effects of Nuclear Weapons. U.S. Government
Printing Office, Washington D.C..
Greenstreet, S., Ngo, H., Gladman, B., January 2012. The orbital distribution of
near-Earth objects inside Earth’s orbit. Icarus 217 (1), 355–366.
Harris, A., 2015. The population of near-Earth asteroids. Icarus 257, 302–312.
Hills, J., Goda, M., 1998. Damage from the impacts of small asteroids. Planet. Space
Sci. 46 (2-3), 219–229.
Hills, J., Goda, M., 1993. The fragmentation of small asteroids in the atmosphere.
Astron. J. 105 (3), 1114–1144.
Ivanov, B., Deniem, D., Neukum, G., 1997. Implementation of dynamic strength models into 2D hydrocodes: applications for atmospheric breakup and impact cratering. Int. J. Impact Eng. 20 (1-5), 411–430.
Johnston, C., Samareh, J., Brandis, A., 2015. Aerothermodynamic characteristics of
16–22 km/s Earth entry. 45th AIAA Thermophysics Conference, AIAA AVIATION
Forum [Online] Available: http://dx.doi.org/10.2514/6.2016-1514.
Kring, D., 1997. Air blast produced by the meteor crater impact event and a reconstruction of the affected environment. Meteorit. Planet. Sci. 32 (4), 517–
530.
Laurence, S., 2006. Proximal Bodies in Hypersonic Flow. California Institute of Technology.
Macke, R., 2010. Survey of Meteorite Physical Properties: Density, Porosity and Magnetic Susceptibility PhD Thesis. Dept. Phys., University of Central Florida, Orlando, FL.
Mainzer, A., et al., 2011. Neowise observations of near-Earth objects: preliminary
results. Astrophys. J. 743 (2).
Melosh, H., 2007. Physical effects of comet and asteroid impacts: beyond the
crater rim. In: Comet/Asteroid Impacts and Human Society. Springer, Berlin,
pp. 211–224.
Morrison, D., Venkatapathy, E., 2017. Asteroid generated tsunami: summary
of NASA/NOAA National Aeronautics and Space Administration workshop
NASA-TM-219463.
Morrison, D., 1992. The Spaceguard Survey: Report of the NASA International
Near-Earth-Object Detection Workshop NASA-TM-107979, NAS 1.15:107979.
Morrison, D., Chapman, C., Slovic, P., 1994. The impact hazard. In: Hazards Due to
Comets and Asteroids. Tuscon: Univ. of Ariz. Press, pp. 59–92.
Motiwala, S., Mathias, D., Mattenberber, C., 2015. An integrated physics-based risk
model for assessing the asteroid threat. Probabilistic Safety Assessment Conference.
NASA Authorization Act of 2005. December 30, 2005, U.S. House. 109th Congress,
1st Session., Public Law 109-155. Available online: https://www.gpo.gov/
fdsys/pkg/PLAW-109publ155/pdf/PLAW-109publ155.pdf.
National Research Council, 2010. Defending Planet Earth. National Academies Press,
Washington, D.C.
Opik, E., 1958. Physics of Meteor Flight in the Atmosphere. USA: Interscience, Mineola, New York.
Passey, Q., Melosh, H., 1980. Effects of atmospheric breakup on crater field formation. Icarus 42, 211–233.
Popova, O., et al., 2013. Chelyabinsk airburst, damage assessment, meteorite recovery, and characterization. Science 342 (6162), 1069–1073.
Popova, O., et al., 2011. Very low strengths of interplanetary meteoroids and small
asteroids. Meteorit. Planet. Sci. 46 (10), 1525–1550.
Register, P., Mathias, D., Wheeler, L., 2017. Asteroid fragmentation approaches for
modeling atmospheric energy deposition. Icarus 284C, 157–166. doi:10.1016/j.
icarus.2016.11.020.
Reinhardt, J., Chen, X., Liu, W., Manchev, P., Pate-Cornell, M., 2016. Asteroid risk assessment: a probabilistic approach. Risk Anal. 36 (2), 244–261.
Rumpf, C., Lewis, H., Atkinson, P., 2016. On the influence of impact effect modelling for global asteroid impact risk distribution. Acta Astronautica 123,
165–170.
Shoemaker, E., 1961. Interpretation of lunar craters. In: Physics and Astronomy of
the Moon. Academic Press, New York, pp. 283–359.
Spohn, T., Breuer, D., Johnson, T., 2014. Encyclopedia of the Solar System, 3rd ed.
Elsevier, Waltham, MA, USA.
Stokes, G., et al., 2003. Study to Determine the Feasibility of Extending the Search
for Near-Earth Objects to Smaller Limiting Diameters. National Aeronautics and
Space Administration.
D.L. Mathias et al. / Icarus 289 (2017) 106–119 119
Toon, O., Zahnle, K., Morrison, D., Turco, R., Covey, C., 1997. Environmental perturbations caused by the impacts of asteroids and comets. Rev. Geophys. 35 (1),
41–78.
Vasilyev, N., 1998. The Tunguska meteorite problem today. Planet. Space Sci. 46 (2),
129–150.
Wheeler, L., Register, P., Mathias, D., 2017. A fragment-cloud approach for modeling
asteroid breakup and atmospheric energy deposition. Icarus, in press.
Yoshinaka, R., Osada, M., Park, H., Sasaki, T., Sasaki, K., 2008. Practical determination
of mechanical design parameters of intact rock considering scale effect. Eng.
Geol. 96 (3-4), 173–186.
</A probabilistic asteroid impact risk model: assessment of sub-300 m
impacts>


<NEOviz: Uncertainty-Driven Visual Analysis of Asteroid Trajectories>
NEOviz: Uncertainty-Driven Visual Analysis of Asteroid Trajectories
Fangfei Lan , Malin Ejdbo , Joachim Moeyens , Bei Wang , Anders Ynnerman , Alexander Bock
Fig. 1: NEOviz showing a prediction of (99942) Apophis’s impact in 2029. The blue Uncertainty Tube is a bundled representation of the
trajectory samples of which only few (green) impact Earth. The orange cross-section shows the extent of all impacting trajectories. The
reference grid indicates a size of 50,000 km. The inset shows the predicted locations for all impacting trajectories on Earth.
Abstract—We introduce NEOviz, an interactive visualization system designed to assist planetary defense experts in the visual analysis
of the movements of near-Earth objects in the Solar System that might prove hazardous to Earth. Asteroids are often discovered using
optical telescopes and their trajectories are calculated from images, resulting in an inherent asymmetric uncertainty in their position
and velocity. Consequently, we typically cannot determine the exact trajectory of an asteroid, and an ensemble of trajectories must
be generated to estimate an asteroid’s movement over time. When propagating these ensembles over decades, it is challenging to
visualize the varying paths and determine their potential impact on Earth, which could cause catastrophic damage. NEOviz equips
experts with the necessary tools to effectively analyze the existing catalog of asteroid observations. In particular, we present a
novel approach for visualizing the 3D uncertainty region through which an asteroid travels, while providing accurate spatial context in
relation to system-critical infrastructure such as Earth, the Moon, and artificial satellites. Furthermore, we use NEOviz to visualize
the divergence of asteroid trajectories, capturing high-variance events in an asteroid’s orbital properties. For potential impactors, we
combine the 3D visualization with an uncertainty-aware impact map to illustrate the potential risks to human populations. NEOviz was
developed with continuous input from members of the planetary defense community through a participatory design process. It is
exemplified in three real-world use cases and evaluated via expert feedback interviews.
Index Terms—Astrovisualization, planetary defense, uncertainty visualization, asteroids.
1 INTRODUCTION
Our Solar System is filled with asteroids and comets, making it inevitable that a large object will impact Earth. An estimated 1 trillion
objects orbit the Sun, out of which around 1.3 million asteroids and
4,000 comets have been discovered so far. The number of newly identified small bodies has been growing exponentially in recent years [11].
Of particular interest are the near-Earth objects (NEOs). A NEO is
a comet or an asteroid whose trajectory brings it close to Earth’s orbit around the Sun. One of the most significant NEOs in history is
the Chicxulub impactor, an asteroid with an approximate diameter of
10 km. It struck Earth about 66 million years ago, ending the age of the
dinosaurs, marking one of the most devastating events in the history
of life on Earth [33]. Currently, we know of roughly 75 asteroids of
comparable or greater size that would have civilization-ending consequences if they were to collide with Earth. More recently, the Tunguska
• Fangfei Lan and Malin Ejdbo contributed equally.
• Fangfei Lan and Bei Wang are with the University of Utah.
E-mail: {fangfei.lan, beiwang}@sci.utah.edu
• Malin Ejdbo, Anders Ynnerman, and Alexander Bock are with Linköping
University.
E-mail: {malin.ejdbo, anders.ynnerman, alexander.bock}@liu.se
• Joachim Moeyens is with the University of Washington.
E-mail: moeyensj@uw.edu
event in 1908 flattened around 2,000 km2
of forest and was caused by
an asteroid estimated to be 50-60 m [10]. In 2013, the 20 m Chelyabinsk
meteor exploded over Chelyabinsk, injuring over 1,600 people.
Planetary defense is concerned with detecting potential impactors
and devising strategies for their mitigation. The risks posed by these
potential impactors are particularly prevalent if the asteroids airburst
in or near densely populated areas, which can lead to catastrophic loss
of life. We develop NEOviz as a visualization system that enables
uncertainty-driven analysis of asteroid trajectories using a participatory design process that involves collaboration with planetary defense
experts from the B612 Foundation. The foundation’s mission is to
discover and track potentially hazardous NEOs, raise awareness of their
risks, and devise strategies to mitigate their impact.
When observing an asteroid, it is much easier to infer its (x, y) position from a telescope image than the depth in the line-of-sight direction
z. This introduces inherent asymmetric uncertainty into the orbit fitting
and propagation process, which is challenging to visualize. The limited
number of observations, combined with the fact that most observations
can only occur at night and can take place over several months, leads
to a scarcity of data when accurately reconstructing the trajectory of an
asteroid, introducing another source of uncertainty. While the quantification of such uncertainty using, for example, covariance matrices, at
specific time points and corresponding mitigation strategies has been
well-researched [30, 39], the visualization of its progression over time
and its impact on Earth remains largely unexplored. An effective and
accessible tool for communicating the complex, time-varying 3D uncertainty could be invaluable not only for scientific research but also
arXiv:2411.02812v1 [astro-ph.EP] 5 Nov 2024
for public outreach and education.
Designing an informative uncertainty representation for asteroid
orbit trajectories presents several challenges, including visual clutter,
time-varying uncertainty, and significant variability in both time and
space. Since the impact probabilities for any given asteroid are fairly
low (for example 2.7% for (99942) Apophis), we can only obtain a representative ensemble of an asteroid’s potential trajectories by drawing a
large number of samples from its uncertainty distribution. This process
creates inevitable visual clutter and heavy computational cost, as many
thousands of orbits must be displayed, with only a few trajectories
impacting Earth. Additionally, as the uncertainty evolves over time, it
is crucial to visualize its progression instead of a collection of static
snapshots, which requires a smooth transition between consecutive
time steps. Finally, there is large variability in both the magnitude of
the uncertainty and the corresponding time period as some uncertainty
propagation could span several decades, while others only involve a few
hours of data. For some asteroids, their orbit uncertainty could exceed
the diameter of the Earth, while for others, their impact point can be
predicted accurately. Therefore, we strive to design our visualization
system to accommodate and adapt to these wide ranges of scenarios.
In response to these challenges, we developed NEOviz, a visualization system that supports the analysis of individual asteroid trajectory
uncertainties. NEOviz is integrated into the astrovisualization engine
OpenSpace [6], which provides the necessary spatial context for asteroid observations. The resulting visualization system offers capabilities
surpassing currently available methods, enabling experts to gain important scientific and operational insights into an asteroid’s potential
movements. The major components of NEOviz include:
• Given an ensemble of an asteroid’s trajectories, we present the
Uncertainty Tube, a spatiotemporal representation of the volume
encompassing all sampled asteroid trajectories. It provides a timevarying visualization of the physical variability of the ensemble.
• At the cross-sections of the Uncertainty Tube, Interactive cutplanes shows the sample distribution at particular time steps. For
each time step, the user can browse the cut-planes and examine
each distribution in detail along the Uncertainty Tube.
• The Impact Map visualizes the likelihood of impact locations on
Earth together with an estimated risk analysis for each location.
The rest of the paper is structured as follows. Sec. 2 discusses prior work
on ensemble and uncertainty visualization focusing on curve-based and
3D uncertainties, and then surveys NEO-related visualization tools.
Sec. 3 provides technical knowledge about the unique characteristics
of this uncertain data. We describe our system components in detail
in Sec. 4 and three case studies in Sec. 5. Finally, we summarize the
results of four qualitative interviews with planetary defense experts in
Sec. 6 and discuss limitations and future work in Sec. 7.
2 RELATED WORK
Uncertainty quantification is a vibrant field of study in astronomy due
to the inherent error margins introduced in the data acquisition process.
Thus far, visualization has, however, only rarely been used to convey
this uncertainty in 3D, because of the data complexity and the large
variability in both time and space. In this section, we first review
relevant prior works on ensemble and uncertainty visualization. We
then discuss existing methods for analyzing the trajectories of NEO.
2.1 Ensemble and Uncertainty Visualization
Significant advances have been made to ensemble visualization due
to the increased availability of data, which is often a result of running
multiple instances of the same physical simulation [36]. Uncertainty visualization is a common technique to effectively communicate the corresponding uncertainty. Several surveys detailing its development and applications are available [17,27–29]. Recently, Wang et al. [36] surveyed
visual analysis techniques on ensemble simulation data and partitioned
scientific visualization techniques of spatial data into four categories:
point-oriented, curve-oriented, surface-oriented, and volume-oriented
approaches. Highly relevant to our work, we next discuss curve-based
uncertainty techniques, especially in the broader context of 3D and
temporal uncertainties; an area which in itself has seen limited attention
from the research community thus far.
Zhang et al. [40] used variable spatial spreading and spaghetti plots
to show the uncertainty in the attribute variable dimension. Their visualization was used to identify features of interest and perform comparative
analysis in a weather forecast ensemble. Inspired by conventional boxplots, Whitaker et al. [38] introduced contour boxplots by computing
order statistics on an ensemble of contours. The contour boxplot was
then further enhanced through the use of curve boxplots [23] as an
extension to streamlines and pathlines. Liu et al. [22] moved away
from traditional error cones and proposed an implicit uncertainty visualization by selectively sampling from the original ensemble and
reconstructing a spatially well-organized ensemble. They applied their
method to tropical cyclone forecast tracks and demonstrated that their
visualization could assist scientists in the exploration and estimation of
storm damages while preventing visual confusion. Analogous to their
approach, we employ a combination of a spatially distributed ensemble,
an extension of an error cone into a time-varying 3D structure, and the
visualization of individual trajectories to more effectively convey the
potential hazards associated with a particular near-Earth object.
Visualizing uncertainties in 3D is known to be challenging.
Guo et al. [14] developed eFLAA to characterize the variations among
3D ensemble flow fields. They incorporated a 3D navigation and a
timeline view to provide an overview and regional comparison capabilities based on the user’s selection. A 3D uncertainty band was also used
to visualize geospatial uncertainties. Fersl et al. [12], on the other hand,
visualized uncertainty in vector field ensembles using confidence lobes,
which are constructed using dimensionality reduction and clustering
on the ensemble. To prevent disastrous underground utility strikes,
Liet al. [21] used multipatch surface models to construct 3D probabilistic uncertainty bands that enclose the true utility lines while accounting
for positional uncertainties. These uncertainty representations inspired
the design of the Uncertainty Tube in NEOviz.
Zhang et al. [41] used streamtubes to visualize 3D diffusion tensor
MRI data, where streamlines are grouped to generate a tube representation. The main differences to our approach are that we must consider
all sampled trajectories due to the asteroids’ small impact probabilities, whereas it was sufficient to characterize the main behavior of the
streamlines in their data by selecting a small subset of trajectories; we
need to perform time stitching and ensure smooth transition during the
time steps while their method shows a static tube visualization.
In astrophysics, few visualizations exist to facilitate exploration of
uncertainty. Bock et al. [7] proposed a system to study the uncertainty
of space weather simulations. They provided a multi-view visualization providing comparative analysis capabilities, the exploration of
time-varying components in the data, and a volume rendering of the
simulations. In a recent survey on the visualization of astrophysics,
Lan et al. [19] identified uncertainty visualization as a significant challenge in the field, offering numerous research opportunities. In this
work, we present a system that represents an effort to address this gap.
2.2 Visualization in NEO Risk Analysis
The analysis of NEO risk is a well-explored subject with considerable
effort dedicated to quantifying the implications of asteroid orbit uncertainty. However, visualization techniques in the field have remained
rudimentary. Most studies use bar charts, line plots, and scatter plots as
the main methods of communicating their scientific insights.
Yeomans et al. [39] were among the first to raise concerns about
(99942) Apophis, discussing keyhole events in both 2029 and 2036,
and considering potential deflection strategies for near-Earth asteroids.
They introduced a compelling “keyhole map”, an adaptation of a line
plot, that is generated by systematically probing the current uncertainty
region of the trajectory and propagating each potential trajectory forward. The map offers an overview of how variations in the current
orbit could result in catastrophic consequences decades later. Paek et al.
then further elaborated on the existence of multiple keyhole events for
Apophis in 2036 and investigated various deflection strategies [25, 26].
Wheeler et al. [37] analyzed which asteroid properties and entry
parameters cause the most damage in case of an impact. They used
tornado plots and histograms to conduct a comparative analysis of
the risk uncertainties for four nominal asteroid sizes ranging from 50–
500 m in diameter and concluded that the impact location contributes
most to impact risk uncertainty, followed by its size and velocity.
Similar to our goal of facilitating asteroid risk analysis, Rumpf
et al. [31] proposed a probabilistic visual representation of impact
corridors. They also incorporated uncertainty information by using a
Monte Carlo method to sample the orbital solution state space based
on the covariance matrix and propagating the orbits forward to compute impact locations. They visualized the spatial impact probability
distribution of each impact location using a 2D scalar field and scatter
plots. In a later work, Rumpf et al. [32] investigated how uncertainty
and risk were affected by a planned deflection of an uncertain asteroid.
In their work, they visualized the impact corridor with a 2D scatter plot.
Norlund et al. [24] provided a more reliable prediction of human casualties and infrastructure risks caused by a NEO impact in NEOMiSS.
They developed a behavior-based evacuation model by combining the
physical effects of a potential impact, the historical insights into the
impact region such as natural hazards and local building properties,
and crowd-sourced information regarding transportation infrastructure.
They overlaid a 2D map with the population density of the impact
region and color-coded each area based on the predicted human casualties. The Impact Map in NEOviz is inspired by this approach. However,
we extend it to first show the location within its 3D context, providing
the ability to show other georeferenced maps simultaneously.
These state-of-the-art visualizations in planetary defense demonstrate the usage of effective visualizations on their own. However,
there is currently a severe lack of a contextualized integration of these
disparate techniques. The ability to investigate the underlying uncertainty of impact locations strongly depends on the ability to inspect
the NEO trajectories in a time-varying 3D space to gain the necessary
insights. NEOviz aims to address this by combining improved versions
of several visualization techniques into a common reference frame and
simultaneously supporting interactive exploration.
3 BACKGROUND
This section describes the domain background focusing on the characteristic nature of NEO orbits and the uncertainty associated with the
data acquisition and orbit propagation. This section also describes how
the input data for NEOviz was created from these initial observations.
While humanity’s ability to discover asteroids has drastically improved, there is still much left to do. The vast majority of asteroid
discoveries have so far been contributed by ground-based, optical surveys. Starting in 2025, the Vera C. Rubin Observatory will begin a
10-year survey called the Legacy Survey of Space and Time [15], which
is estimated to contribute over 5 million new asteroid discoveries, nearly
quadrupling the known asteroid population. This increase in data exacerbates the need for a visual analytical tool to systematically examine
the new discoveries [15, 16]. Additionally, space-based telescopes not
only serve to characterize the sizes of asteroids but can also discover
asteroids that are not visible from the ground.
Most asteroids and comets are located in the Main Asteroid Belt
between Mars and Jupiter or the Kuiper belt beyond Neptune. Over
short time spans, typically a few weeks, their motion on perturbed
Keplerian orbits can be approximated with simple 2-body dynamics
due to the Sun’s overwhelming gravitational influence. However, over
longer time spans, the gravitational influence of the planets and other
massive objects needs to be accounted for to accurately predict the
positions of asteroids. In addition to gravitational effects, such as
collisions in the Main Belt or close encounters with a planet, asteroids
and comets may also be subject to non-gravitational forces such as
outgassing, the Yarkovsky effect [35] that introduce high levels of
uncertainty and make long-term predictions challenging.
Asteroids are typically discovered at the magnitude limit of astronomical surveys and near opposition where they appear brightest. The
orbital uncertainty within the arc of observations is well-approximated
by a 3D ellipsoid, where the uncertainty of the orbit derived from observations is commonly represented as a 6D covariance matrix. The
six dimensions correspond to the orbital parameters of an asteroid,
namely its Cartesian coordinates (x, y,z) and the velocity in each direction (vx, vy, vz). The diagonal of the covariance matrix depicts the
variance of each orbital parameter. Beyond the arc of observations, the
uncertainty region instead becomes extended or elongated (referred to
as a “banana-oid”). This elongation occurs due to Keplerian shear, in
which orbits closer to the Sun will possess greater velocities than orbits
further from the Sun. Propagated over time, the variation in distance
from the Sun, and by extension, the variation in velocity, causes the
uncertainty region to become elongated or extended.
Further observations through observation campaigns or via precovery [5] in existing archives refine the orbital parameters of the object
and influence its uncertainty. Impact probability studies typically involve sampling the covariance matrix using Monte Carlo methods for
variant orbits and propagating these forward. The ratio between impacting orbits over the total number of variants is then reported as the
object’s impact probability. For this sampling technique, since each
variant is individually propagated, Monte Carlo sampling is both the
most accurate but also the most computationally expensive to perform.
A global effort, coordinated in part by the United Nations Office
for Outer Space Affairs, NASA, ESA, JAXA, and others, aims to identify and mitigate any potential future asteroid impacts. The Asteroid
Institute, a program of the B612 Foundation, is a research institute
that participates in the planetary defense effort by building the tools to
process the next generation of astronomical data. One of these tools
is the Asteroid Discovery Analysis and Mapping platform (ADAM),
which is designed to enable compute-intensive research in astronomy
and currently has services to perform precovery searches for asteroid
observations, asteroid discovery searches in telescopic data, and Monte
Carlo impact probability calculations. Critical to the planetary defense
effort is the ability to communicate the nature of the impact hazard
to both the planetary defense community and to the general public.
The visualization of orbital uncertainty regions and risk corridors are
necessary tools to facilitate such communication.
For each of the asteroids detailed in the case studies in Sec. 5, we
used mpcq1
to gather observation and submission histories. mpcq allows the user to query the Small Bodies Node replica of the Minor
Planet Center (MPC) database for observations of known asteroids and
details about their submission history. As observatories conduct their
campaigns, bundles of observations are submitted to the MPC typically
at the end of a night of observation. The timestamp associated with each
submission is often defined as the timestamp of the last observation
within the submission. Once a submission of observations is accepted
by the MPC, the orbits of any known objects are updated.
4 SYSTEM
This section outlines the design for our orbit uncertainty visualization
system. We begin by detailing the Uncertainty Tube generation process,
followed by a discussion of the graphical rendering techniques. To
facilitate an efficient scientific discovery process, we also implement
several interactivity features. Lastly, for potential impactors, we create
an Impact Map that displays the predicted impact locations along with
associated risk probabilities.
To determine an orbit from a submission of observations, we use
the software find_orb2
, which is extensively used by the planet defense community [20]. Given a set of observations from a submission,
find_orb determines the best-fit orbit, whose uncertainty is calculated
as a covariance matrix during the least-squares orbit correction optimization. The orbit fitting process for each object is run iteratively
over time, with every new submission including the observations from
previous submissions. The outcomes of the orbit determination process
include the observations, submission history, best-fit orbit derived from
these observations, and the time-varying uncertainty.
In addition to find_orb, we use the adam_core3
library for various
orbit-related computations. adam_core contains the shared set of
utilities used by the services hosted on the ADAM platform and also
1https://github.com/B612-Asteroid-Institute/mpcq
2https://www.projectpluto.com/find_orb.htm
3https://github.com/B612-Asteroid-Institute/adam_core
amateur
observations
submissions
covariance matrix
solar system bodies
satellites
orbit
propagation
impact map
NEOviz
uncertainty tube
professional
observations
Fig. 2: A systematic overview of NEOviz. NEO observations and their associated uncertainties are used to propagate trajectories into the future.
These trajectories are then visualized as a collection using the Uncertainty Tubes and an Impact Map in the case of an impact. The results are
presented in context to other solar system bodies and artificial satellites.
provides access to orbital integrators such as PYOORB [13]. In NEOviz,
we employ adam_core to propagate orbits forward in time and perform
Monte-Carlo sampling using the covariance matrix.
4.1 Uncertainty Tube Generation
For each NEO, the submission data consists of a batch of timestamped
observations. For each submission time t, we obtain a covariance
matrix that represents the uncertainty of the orbit from orbit fitting.
The uncertainty is time-varying and asymmetrical, which means that
the uncertainty values can differ greatly for each direction over time.
The goal is to encapsulate the uncertainty using ellipse slices at each
submission time that are oriented perpendicular to the mean velocity
of the samples. Connecting all ellipse slices results in a 3D tube and
enables the expert to visually inspect it within its spatial context. To
construct the Uncertainty Tube, we introduce a three-step process:
1. Orbit propagation and variant sampling: Compute the best-fit
orbit using the submissions and their associated uncertainties and
sample orbit variants from the resulting covariance matrix.
2. Ellipse slice computation: For each time step, compute an ellipse slice perpendicular to the mean velocity of the samples to
approximate the uncertainty region of the orbit.
3. Time stitching: Find correspondences between neighboring ellipse
slices by sampling a fixed number of points on the perimeter of
each ellipse in a deterministic order. Connect the ellipse slices
into a 3D tube.
An alternative approach to using ellipses could be to use convex hulls,
which can provide a more accurate approximation of the uncertainty
boundary. However, establishing correspondences between convex
hulls during time stitching presents significant challenges, and the
resulting 3D structure during interpolation might not be meaningful
and be more difficult to interpret. We also considered constructing the
ellipses through a 2D Eigendecomposition of the covariance matrix,
but such ellipses may only encompass the majority of the trajectories,
which is insufficient for our purpose.
Orbit propagation and variants sampling. Given the entire submission history of an object, we utilize only the observations available
up to a given time. When an asteroid is first discovered with only a
few observations in the submission history, find_orb fails to identify
an orbit. Therefore, we must determine the initial submission time at
which a sufficient number of observations are available for find_orb to
identify an orbit. Proceeding forward from this point in time, we can
freely select a start and end time for the orbit propagation. For each
subsequent submission time, we find a new best-fit orbit, including the
covariance matrix, considering all past observations.
We introduce two uncertainty representations for the 3D tube: sectioned uncertainty and historical uncertainty. Sectioned uncertainty
applies to events where all submission data is present, whereas historical uncertainty is designed for predictions of asteroid trajectories. For
the sectioned uncertainty representation, as illustrated in Fig. 3 (A), we
start at a submission sub0. We then propagate the associated best-fit
orbit forward to shortly before the next submission time sub1. We call
all such orbits the propagated best-fit orbit. We enforce a time gap (tg
in Fig. 3) between the end of sub0 orbit and the start of the sub1 orbit
to ensure that the changes in size and direction of the tube between
submissions can be clearly seen without introducing instantaneous
changes that would draw unwanted attention. When reaching the last
submission, we propagate the orbit forward to a user-specified end time,
at which point the covariance matrix for this propagated best-fit orbit
is reconstructed using a Monte Carlo sampling method (see Sec. 3).
We then uniformly sample a large number of orbit variants from this
distribution in the same time frame, i.e. from submission times sub0 to
shortly before sub1. Using this procedure, we produce a large sample
of possible orbits that represent the orbit uncertainty at each submission
time until we reach the next batch of observations. The Uncertainty
Tube would be sectioned by the submissions, where the shape of the
tube is recomputed and corrected at each new submission time.
For the historical uncertainty representation, we also start at submission time sub0. However, instead of correcting the propagation with
each subsequent submission, we consider only the observations up to
sub0 for the orbit propagation. As illustrated in Fig. 3 (B), we propagate
forward from sub0 to a user-specified end time without considering
new submissions. Using this orbit and covariance matrix, we repeat the
same procedure as in the sectioned uncertainty scenario and uniformly
sample a large number of orbit variants. At sub1, we repeat the process
sub1 sub0 sub2
sub0
sub1
Sectioned Uncertainty
Historical Uncertainty
A
B
tg tg
sub1  tg sub2  tg
Fig. 3: Illustrations of the two different Uncertainty Tube representations:
(A) sectioned uncertainty and (B) historical uncertainty. The timeline
below each image shows the time ranges for each Uncertainty Tube.
mean velocity
plane
mproj (3d)
n mproj (2d)
A B C D F G
mrot(2d)
mproj (2d) c
E
nss
c0
c1
c2
m0
m1
m2
v0
v1
v2
p0
p1
p2
p
Ellipse Slice Computation Time Stitching
Fig. 4: Ellipse slice computation and time stitching pipeline. At each time step, the trajectories are projected into a 2D space to calculate a minimum
enclosing ellipse, which is used as an uncertainty representation. A consistent reference frame that is co-rotating with the NEO orbits is then used to
connect adjacent ellipses and generate the Uncertainty Tube.
and use the new best-fit orbit to obtain a new set of orbit variants and
create another Uncertainty Tube. The historical uncertainty representation recreates historical events to investigate the predictions made using
only the observations available at a particular point in time.
Ellipse slice computation. To capture the orbit’s uncertainty structure in 3D, we generate ellipse slices that reflect the shape of the orbital
distribution at various time steps. We use an adaptive time sampling
strategy based on the time interval between two consecutive submissions that produces one ellipse slice per submission or at least one slice
per day. A clear change in the Uncertainty Tube at each new submission
is indicated by the time gap, denoted as tg in Fig. 3, defined as a small
percentage of the corresponding submission interval. For the Sectioned
Uncertainty Tube illustrated in Fig. 3 (A), since the submission interval
between sub0 and the end of sub0, namely sub1 −tg, is larger than one
day, we insert an additional ellipse slice in between the starting and
ending ellipse slices. In the case of the Historical Uncertainty Tube
(B), ellipse slices continue to be added as time progresses forward at
an interval of approximately one day. In particular instances where an
expert seeks to closely inspect a certain time period, the temporal resolution is further increased to produce a more refined Uncertainty Tube
section, following the “Overview first, details on demand” mantra [34].
From this process, we obtain a predefined array of time steps.
Fig. 4 (left) provides an overview of our ellipse slice computation
pipeline. At a given time step, the positions of all sampled trajectories
form a point cloud, illustrated in blue in Fig. 4 (A). To represent the
uncertainty, the aim is to compute a 2D ellipse slice for this 3D point
cloud, reducing the problem to a 2D task of finding the minimum
enclosing ellipse of the point cloud in a plane.
We define the center c of the point cloud to be the center of the 3D
minimum volume enclosing ellipsoid. Using the point cloud’s center
c, marked in red in Fig. 4 (A), and the mean velocity direction of the
variants as the normal vector n, we determine a 2D mean velocity plane.
All points in the 3D point cloud are then projected onto the plane in
the direction of their instantaneous velocity vector. In Fig. 4 (A), the
blue points are the original 3D point cloud and the projected points are
drawn in green. These projected points are then transformed from the
mean velocity plane into the x-y plane shown in Fig. 4 (B), where it
is now a 2D problem of finding the minimum area enclosing ellipse
for the projected 2D point cloud. The representative ellipse is found
using the algorithm by Bowman and Heath [8]. Fig. 4 (C) shows an
example ellipse of the 2D point cloud in black. The resulting ellipse,
together with the distribution of the 2D point cloud, is transformed into
a cut-plane at the corresponding time point on the Uncertainty Tube.
Time stitching. We compute an ellipse slice for each of the predefined
time steps. We connect ellipses from adjacent time steps by establishing
correspondences between them. The goal is to uniformly sample the
same number of points from the perimeter of each ellipse and match
the points between any two adjacent slices. This way we prevent the
Uncertainty Tube from twisting during interpolation over time.
In order to establish correspondences, we first need to determine
the position on each ellipse to start the sampling process. We need to
make sure that the starting point is approximately the same for each
ellipse with respect to the center of the ellipse and the location of the
sun. We rely on two vectors that are relatively constant with respect to
the asteroid, the normal of the fundamental plane of the Solar System
nss and the direction from the center of the 3D point cloud ci
to the Sun,
denoted vi for each time step ti
, see Fig. 4 (G). Both the normal and the
directional vectors are computed using NASA’s SPICE library [1, 2].
We take the cross product of these two vectors vi and nss and obtain the
vector mi
. However, mi
is often not on the mean velocity plane. We thus
perform an orthogonal projection of m onto the mean velocity plane,
denoted mpro j(3d) in Fig. 4 (A). As we transform the points to 2D, we
also transform the vector onto the 2D plane, shown as mpro j(2d) in
Fig. 4 (B). We rotate the ellipse for the major and minor axes to align
with the x and y axes and the rotated mpro j(2d), presented as mrot(2d),
is drawn in pink in Fig. 4 (C) and (D). For each ellipse, the sampling
starts at intersection point p of mrot(2d) and the ellipse and uniformly
continues for a fixed number of points in a clockwise direction to result
in the sampled points on the boundary of the ellipse. Fig. 4 (E) shows
the sampled points (orange) in relation to the 3D point cloud (blue).
Finally, we interpolate between corresponding points from adjacent
ellipses, e.g. between p0 and p1, p1 and p2 in Fig. 4 (G).
We encode further the density distribution of the orbit variants as the
surface color of the tube in the Uncertainty Tube rendering stage. This
technique was desired by the domain expert to be able to distinguish
distribution changes over long time periods at a quick glance and to see
whether they are equally distributed or skewed (see Fig. 1, left). For
an ellipse slice, we compute the density information by constructing
a ray from the center of the ellipse passing through each point, and
computing the intersection between the ray and the ellipse. Then, for
each sampled point on the ellipse, we count the number of intersection
points surrounding the sampled point within a given radius r. We assign
this value as the density value of each sampled point on the ellipse. In
Fig. 4 (F), we color each sampled point based on its associated density
value, with the darker colors indicating a higher density. To select the
appropriate r, we take the minimum of two quantities: the length of the
minor axis of the ellipse, and the distance between two adjacent sample
points on the ellipse’s boundary.
4.2 Uncertainty Tube Rendering
Using the ellipse data from the previous step, a tube is generated by
connecting adjacent slices in a triangle mesh to form a 3D structure that
visualizes the known uncertainty of the sampled trajectories (see Fig. 5).
The parameter and the transfer function used to colorize the surface
are controlled by the user. Furthermore, the user can choose to show
the ellipse images created in the previous step on the “lid” of the tube
(Fig. 1 (left)). Lastly, the tube can be rendered as a wire-frame which is
used to inspect individual trajectories that lie on the inside of the tube,
see Fig. 5 for an example.
The user can change the simulation time in NEOviz to simultaneously
control the position of all objects. The extent of the Uncertainty Tube
is only shown until the selected simulation time. Since this selected
time might lie in between two ellipse slices, to display reasonable
information at this time, both the tube’s vertices and the cut-plane
information are interpolated between the preceding and succeeding
time steps using linear interpolation. As the relative time distance
between ellipse slices is fairly small, a linear interpolation was deemed
to be sufficient and not to introduce additional errors.
In the case of visualizing the historical uncertainty, in which multiple tubes are calculated, the domain expert chooses which Uncertainty
Tube to show by providing individualized stop times for each representation. This provides the ability to comparatively study the uncertainty’s
temporal evolution when a specific new submission is added, thus
validating its impact on the orbital characteristics.
4.3 Impact Map
In addition to the trajectories of the NEOs, it is also vital to show
an accurate representation of the location where a potential impactor
might collide with the Earth. As these events are rare, it is beneficial to
provide the geospatial context to the impact locations in relation to the
3D trajectory visualization described above. Traditionally, these maps
are represented as impact corridors which show a trail on Earth where
a NEO might impact. However, these impact corridors rarely display
uncertainties and never include information about the original trajectory
of the sample variant that resulted in a specific impact location.
The trajectories calculated in the previous step are also used to
calculate the predicted impact location on the surface. After all impact
points are collected, these are drawn as circles onto an image using
the equirectangular map projection. The radius of the circle is fixed
for each impact, and is adjustable to represent the potential impact
energy, if such information is available for a given impactor. The
circle furthermore uses a Gaussian falloff to convey the remaining
uncertainty of these impact locations due to atmospheric effects that
are not accounted for in the trajectory calculation. All impact circles
are composited using additive blending to result in a final map that
represents a higher impact likelihood with a higher numerical value.
The impact probability map (see Fig. 1, right) shows areas with a
higher impact probability with a brighter color than areas with individual dots. This indicates where on Earth it would be more or less
likely for the asteroid to impact, giving an approximate of the impact
probability for that location. This map, however, only shows the impact probability irrespective of what is located underneath the impact
location. To visualize the risk to human lives, a second map is generated that highlights areas with large human populations. This map is
generated by combining the impact probability map with the location
of human-made artificial lighting on Earth, which has been shown to
accurately represent the human footprint on Earth [18]. Similar to the
Fig. 5: Toggling the Uncertainty Tube into a wireframe mode enables
the inspection of trajectory variants inside the tube in the case of nested
uncertainty. The blue tube shows the uncertainty of all trajectories, the
orange tube shows ones intersecting Earth.
impact probability, this map is colorized using a transfer function that
maximizes the data’s visibility when projected onto the globe in the 3D
visualization and when composited with existing satellite images.
5 CASE STUDIES
The following three case studies demonstrate various use cases of
NEOviz to facilitate the scientific discovery process of planetary defense
experts. We investigate three asteroids of interest: (367943) Duende,
(99942) Apophis, and 2023 CX1. For each asteroid, we generate 10,000
variant orbits for each best-fit orbit to construct our Uncertainty Tube.
These three case studies are also shown in greater detail in the supplemental material accompanying this paper. The domain expert selected
these asteroids to represent a wide range of challenging scenarios that
the visualizations of NEOviz could address:
Close approach (Sec. 5.1) A close approach with Earth drastically
alters a NEO’s orbital parameters, increasing the uncertainty in unpredictable ways as the resulting changes depend on the distance between
the object and Earth.
Predicted impactor (Sec. 5.2) Predicting the location many decades
into the future is challenging as its covariance matrix quickly becomes
unbounded without further observations. This use case illustrates
NEOviz’s ability to inspect the time-evolution of impactors in 3D, investigate the uncertainty in location over 25 years of orbits, and display
an Impact Map for the small likelihood of an impact.
Imminent impactor (Sec. 5.3) Many NEOs are only discovered hours
or days before their impact. This use case demonstrates NEOviz’s
ability to visualize an object with such a small number of observations
and to aid the inspection of historical uncertainty.
5.1 Close Approach of (367943) Duende
Asteroid (367943) Duende is a small near-Earth asteroid that made
a record-close approach on February 15th, 2013. The flyby was so
close that it traversed inside the geosynchronous satellites, reaching the
closest distance of about 27,700 kilometers from Earth’s surface. Using
NEOviz , we visualize the changes in the uncertainty of orbit prediction
as the number of observations increases and show the progression of
the close approach in relation to the associated uncertainty.
Fig. 7 shows the recreation of the asteroid orbit predictions at two
different times. We use the historical uncertainty representation of
the tube and focus on two submissions, one on March 5th, 2012, one
month after its discovery, and one on February 10th, 2013, five days
before the close approach. The orbit trajectories are propagated forward
from each of these submissions to an end time after the close approach.
Two nested Uncertainty Tubes are generated. Fig. 7 (A) shows both
tubes shortly before the close approach. The larger blue tube represents
the trajectory uncertainty predicted using data from March 5th, 2012.
The uncertainty region is much larger than Earth with a relatively
high probability of intersection. The smaller orange tube represents
an uncertainty region, computed using more recent observations, that
would approach but not intersect Earth. We also observe a more circular
cut-plane and more uniformly distributed points due to a decrease in
uncertainty, indicating a more confident estimation of the asteroid
trajectory. NEOviz effectively showcases the substantial reduction in
uncertainty as the number of observations increases.
To show the progression of the close approach, we generate a sectioned uncertainty tube for the asteroid using all known observations.
We compute a positional uncertainty magnitude at a given time, using
the covariance matrix. We consider the first quadrant of the covariance
matrix, which only represents the positional variance. The three dimensions correspond to the Cartesian coordinates of the asteroid (x, y,z).
The diagonal of the quadrant, (σx,σy,σz) represents the variance in
each of the dimensions, respectively. We compute the positional uncertainty magnitude as the geometric length of that diagonal. This
uncertainty magnitude approximates the radius of the uncertainty region which is correlated to the longest axis of the Uncertainty Tube.
Fig. 6 (A) is a dual-axis plot displaying the inverse correlation between
the uncertainty magnitude in blue and the number of observations in
orange. The three timestamps before, during, and after the close approach are marked in orange on the uncertainty magnitude, which is
(a) Uncertainty magnitude vs number of observations (b) Before close approach (c) During close approach (d) After close approach
Fig. 6: During asteroid (367943) Duende’s close approach on Feb 15th, 2013, the positional uncertainty increases due to the gravitational pull of
Earth, while a second source of uncertainty, a lack of observations, decreases. NEOviz shows the transition accurately as the uncertainty shows a
skew towards the Earth before the closest approach and a more uniform distribution after the encounter. The orange points in (a) correspond to the
times at which figures (b), (c), and (d) were created. The trails around Earth show the Geostationary satellites at a distance of 36,000 km.
Fig. 7: Nested Uncertainty Tubes for asteroid (367943) Duende as
propagated from submissions on March 5th, 2012 (blue), and February
10th, 2013 (orange), respectively. (B) is a zoom-in view of (A). The
distance between grid cells is 7,500 km.
also shown in Fig. 6 (B), (C), and (D), respectively, where the yellow
trails show the geostationary satellites. The uncertainty is large before
the close approach due to a lack of observations. We observe an interesting phenomenon on the cut-plane of the Uncertainty Tube, where
most trajectories accumulate on the side of the tube in the direction
of the Earth, indicating the asteroid’s movement approaching Earth.
The uncertainty decreases as the asteroid moves closer to Earth. The
Uncertainty Tube reduces in size and the cut-planes circularize as the
distance uncertainty decreases drastically with increased numbers of
observations. This perception corresponds exactly to the findings in the
uncertainty plot. The Uncertainty Tube enables a visual inspection of
the changes in asymmetric uncertainty in more detail.
5.2 Historical High Probability Impact of (99942) Apophis
For this case study, we recreate a historical prediction of the high
probability impact event of the asteroid (99942) Apophis, estimated to
be 300 to 400 meters along its longest axis. Upon its discovery in June
2004, it quickly gained notoriety as it was recognized as potentially
hazardous based on the observations of December 20th, 2004, where
Fig. 8: A point cloud visualization showing 1000 possible trajectories
for Apophis, moments before impact on April 13th, 2029. The impactor
trajectories, marked in red, are visibly more affected by Earth’s gravity.
JPL predicted the asteroid’s potential impact on Friday, April 13th, 2029
with an impact probability of 0.02%. However, the impact probability
continued to increase, peaking on December 27th at a 2.7% likelihood
of impact in 2029 [9]. Although the impact in 2029 was later ruled out
following further observations, it remains a significant historical event
with potentially daunting implications. We recreate this high-impact
probability prediction using NEOviz to show the possible trajectories of
Apophis given the available data on December 27th, 2004 (see Fig. 8).
We compute a historical uncertainty tube for Apophis from the
last submission on December 27th, 2004 to shortly after the predicted
impact time on April 13th, 2029. We sample 10,000 trails from the
uncertainty distribution at the initial submission and propagate forward.
We observe a slight decrease in uncertainty hours before the predicted
impact time and an immediate and drastic increase in uncertainty after approaching Earth. Shown both in the uncertainty magnitude plot
in Fig. 9 (A) and the NEOviz visualization in Fig. 10. The Uncertainty
Tube decreased in diameter from (A) to (B) but increased in size significantly from (B) to (C). Out of these 10,000 trails, we identify 115 trails
that impact Earth, indicating an impact probability of 1.15%, comparable to the historical record [9]. We generate a second Uncertainty Tube
using only the impact trails. We observe in Fig. 11 (top) that at the
same point in time, the impactor tube is further ahead while the larger
tube lags behind. This indicates an interesting phenomenon that the
impactor trajectories traveled further on average than the non-impactor
trajectories. As we move forward in time to the predicted impact time
around 21:00 UTC on April 13th, 2029, we show both a nested tube and
the impactor tube alone, see Fig. 11. We witness a small elevation in
the sub-tube of the impactor trajectories, representing the perturbation
in the asteroid orbit caused by Earth’s gravitational pull.
We further investigate the cause of this phenomenon. Since the point
distribution on the cut-planes is reflected in the surface color of the
Uncertainty Tube, we can change the color map on the surface to reveal
the changes in the point distribution more clearly over time, as shown
(a) Uncertainty magnitude (b) Impactor vs non-impactor samples
Fig. 9: The uncertainty magnitude for Apophis drastically increases after
the predicted impact in 2029. The cut-plane shows the cluster of samples
that will impact (orange) and miss (blue) Earth.
(a) Before close approach (b) During close approach (c) After close approach
Fig. 10: (99942) Apophis impact event on April 13th, 2029 at time 21:43:59 UTC. Note that the data used to generate this image comes from
observational data in 2004 when the impact probability for the object was at its highest. This does not reflect what will happen in 2029 as additional
observations since 2004 have ruled out an impact. Due to the gravitational influence of Earth, the distribution of trajectory samples inside the
Uncertainty Tube shifts as it approaches Earth. The close encounter causes a rapid increase of the uncertainty (see also Fig. 9, A). The grid lines of
the gray reference grid have a size of 30,000 km. The yellow lines around Earth represent the geostationary satellites at a distance of 36,000 km.
in the large blue tube in Fig. 1. We observe that the surface colors
change drastically approximately 2 days before the impact date. At the
time of the predicted impact, the points on the cut-plane diverged from
an elliptical Gaussian to the plot shown in Fig. 9 (B). We see that the
impact trails (orange) are clustered together in relation to all sampled
trails (blue). This change in shape of the impactor tube suggests that
the asteroid is near Earth for some time before the impact, allowing the
Earth’s gravity to have a large effect on its orbit. We can also deduce
that the asteroid has low impact velocity or v∞ before it was accelerated
by Earth’s gravity. Our finding is confirmed by JPL’s Small-Body
Database4
, where the v∞ for Apophis on the impact date is predicted to
be 5.84 compared to an average value of about 20. We show this effect
more directly using a time-varying point cloud visualization of 1,000
trajectories, see Fig. 8. The red points are the impactor trajectories and
the green are regular trajectories. We see that the impactor trajectories
are pulled towards Earth and diverge from the rest of the ensembles.
4https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/
Fig. 11: The top image combines two Uncertainty Tubes to show a
combination of all trajectories (blue) and only impactors (orange). The
Uncertainty Tube is wider than Earth as impacting trajectories can lead
or lag behind the mean velocity plane and thus arrive at Earth at different
times. The bottom image displays the impactor tube only.
For all 115 potential impact trajectories for Apophis, we compute
their impact locations and generate an Impact Map (see Sec. 4). We
present the impact probability (Fig. 1) and the population risk views
(Fig. 12). In the impact probability view, we see many potential impact
points in Africa, a few in the Indian Ocean, and Europe. The brighter
impact points, reflecting higher impact probability, are mostly located
close to Madagascar and Tanzania. In the population risk view, the
brightness of the impact points is determined by the density of the
population at the impact locations. We see the brightest points in
Tripoli, Patna, Riyadh, Tashkent, and Bangkok, corresponding to areas
of high population. This view offers compelling visual cues of urgency
that could potentially aid in the development of evacuation strategies to
mitigate disastrous outcomes of an asteroid impact.
5.3 Imminent Impactor 2023 CX1
2023 CX1 was a small asteroid discovered at the Konkoly Observatory
in Hungary on February 12th, 2023 around 21:20 UTC that entered
Earth’s atmosphere over northwestern France six hours later at approximately 2:59 UTC. The asteroid created a bright fireball across the
English Channel, which was visible in France, the United Kingdom,
Germany, and Belgium [3]. A previous visualization of the impact
corridor provided by the European Space Agency displays a short band
moving towards Normandy, France5
. NEOviz shows a comparable
result, augmented with additional information and conveyed with enhanced visual impact. Due to the extremely short window between
initial discovery and impact, only a few hundred observations are available for 2023 CX1. This results in a highly asymmetrical historical
5https://neo.ssa.esa.int/past-impactors/2023cx1
Fig. 12: Risk assessment for (99942) Apophis in 2029 showing the
convolution between the likelihood of an impact with population density.
uncertainty tube shown in Fig. 13. Using all available knowledge, we
show the tube from 02:38 to 03:40 on Feb 13th, 2023. As the asteroid
is an imminent impactor, the impact location is constrained to a small
region as the uncertainty does not have time to disperse over a larger
surface area before impacting. Fig. 13 shows the Uncertainty Tube and
the Impact Map created using NEOviz. The Uncertainty Tube has a
narrow shape, representing a relatively high uncertainty in the depth
direction. Our predicted impact region is close to the coast in Normandy, France, aligning with the other predictions. NEOviz presents
an interactive tool to simultaneously display the asteroid trajectory and
the location of the imminent impactor.
6 EXPERT FEEDBACK
The development of NEOviz followed a participatory design process
in collaboration with a domain expert from the B612 foundation. To
more systematically assess the effectiveness of our visualization, we
conducted a series of in-depth interviews with five experts, who had
no prior experience with NEOviz. These experts were co-founder and
director of the B612 Foundation, Dr. Ed Lu, and renowned researchers
in asteroids and small solar system bodies, Dr. Aren Heinze, Dr. Mario
Juric, and Dr. Siegried Eggl. All experts have extensive expertise on ´
NEOs and their associated uncertainty, as well as substantial experience
using visualization in their research; one also has experience developing
new visualization tools. During the interviews, we presented a tutorial
of NEOviz, followed by an interactive demonstration of the Apophis
case study. The experts then provided feedback on the effectiveness of
each component of the visualization tool and offered suggestions for
improvement. We refer to the experts as E1-E5, in no particular order.
We received positive feedback and encouragement for NEOviz, especially for its novel approach to visualizing uncertainty. E1 acknowledged that we are tackling a “tough visualization task”, and remarked,
“I haven’t seen anything like this, and I think it is great”. Both E3 and E4
discussed the difficulties they face in visualizing trajectory uncertainty,
noting that NEOviz is a unique tool with the potential to “be a game
changer” for the planetary defense community.
For the Uncertainty Tube and the cut-planes, E1 commented that they
are superior to the point cloud visualization as one can clearly see the
shape of the complex 3D structure in a time-varying setting. E3 viewed
the cut-planes as an effective way to “provide a 2D snapshot of a fairly
complex 3D structure” and appreciated the changing of tube shape as
Fig. 13: The imminent impactor 2023 CX1 on 2023 February 13 02:59:21
UTC is visualized with the Uncertainty Tube. The top image shows an
orange tube that is highly asymmetrical. The lower image displays the
impact location. As the impact occurred during the night, Earth was
artificially lit to make the surrounding geography visible.
new observations arrived. E2 and E4 emphasized that the distribution
of the points is useful for further analysis. E2 found the elevation
in the tube before the Apophis potential impact “fascinating”, saying
that “we need visualization like this to see” such features. However,
experts also provided valuable suggestions for improvement. E1 and E2
expressed reservations about removing the temporal aspects from the
visualization, as the cut-planes depict the projection of the trajectories
instead of their true locations. E2 proposed adding an option for a cutsurface that respects the true time of each trajectory, in addition to the
cut-planes. E4 suggested providing a statistical view on the cut-planes
by showing a kernel density estimation map of the point distribution.
Experts responded positively towards the impact map, with E4 stating it was an excellent idea to convolute population density on top of
the impact locations. E1 suggested offering users more flexibility in
adjusting the size, color, and brightness of the impact markers. E2 and
E4 expressed interest in using an upsampling approach to form more
defined impact corridors instead of scattered impact points. However,
this is currently beyond our computational and rendering capabilities.
NEOviz generated several surprising new insights and inspired many
promising research directions. E2 and E3 were both fascinated by the
elevation in the impactor tube, and expressed curiosity about the extent
of trajectory scrambling over the 15-year propagation period. E3 found
the point cloud visualization particularly engaging, especially where
the impactor orbits are visibly pulled towards Earth, stating “that’s
something you don’t often think about, and you obviously see it here.”
E3 also highlighted potential applications of NEOviz in visualizing
spacecraft trajectories, and the possibility of generalizing it to star
orbits, while cautioning about associated complexities. E4 discussed
applying NEOviz to space collision problems, such as space debris and
collisions within the asteroid belt. On a broader note, E2 appreciated
NEOviz as a tool for unscheduled exploration of asteroid data, stating,
“I spend a lot of time staring at data and just exploring it. Sometimes it
feels like a waste, but I think it is actually really productive in the long
run. And this seems like a tool one can use for that kind of exploration.”
7 CONCLUSION AND FUTURE WORK
In this paper, we present NEOviz, a system enabling planetary defense
experts to analyze the time-varying uncertainty in the trajectories of
NEOs. Our approach introduces two key visualizations: the Uncertainty Tube and the Impact Map. The Uncertainty Tube provides a novel
time-varying 3D uncertainty representation of the asteroid trajectories,
allowing users to explore the temporal evolution of the orbital uncertainty in context with other objects in the Solar System, such as Earth
and satellites. The Impact Map, on the other hand, displays potential
impact locations, combined with population risk, on a 3D model of
Earth using the GlobeBrowsing system [4]. By integrating these visualizations into a common reference frame, NEOviz enables planetary
defense experts to analyze the uncertainty evolution of NEO trajectories
to an extent that was unattainable with previous techniques. NEOviz addresses unique challenges in working with asteroid data, which exhibits
large variability in both time and space due to the diversity in asteroid
scenarios. We demonstrate the capabilities of NEOviz with three asteroids: (367943) Duende, (99942) Apophis, and 2023 CX1, and further
evaluate our system with expert feedback interviews.
While the expert interviews showed that NEOviz has significant
promise as a tool for analyzing asteroid trajectories, it has some limitations. Projecting the 3D point cloud onto a 2D plane results in
information loss and inaccuracies in the uncertainty representation.
Additionally, the current system does not support displaying the numerical data about the uncertainty tube or the asteroid. Addressing these
limitations could be valuable directions for future work. Moreover, the
expert interviews highlighted the scientific benefits of providing statistical views of tbe uncertainty representation by using kernel density
estimation. The orbit variant generation process could also be further
improved.
In light of new telescopes providing an abundance of asteroid data,
we would like to extend the system to visualize multiple asteroids simultaneously, and to enable comparative analysis of NEOs on a population
level rather than individually. Furthermore, we see the potential of
using NEOviz for public outreach. The system could be used to create
immersive experiences in planetariums, raise awareness of NEOs, and
ultimately garner public support for funding future missions aimed at
safeguarding Earth from catastrophic impacts.
ACKNOWLEDGMENTS
This work was supported by the Knut and Alice Wallenberg Foundation
under grant KAW 2019.0024 and through the WISDOME project,
NASA under grant NNX16AB93, the Swedish e-Science Research
Centre, the DIRAC Institute in the Department of Astronomy at the
University of Washington, and the Asteroid Institute, a program of
B612. The DIRAC Institute is supported through generous gifts from
the Charles and Lisa Simonyi Fund for Arts and Sciences, and the
Washington Research Foundation. The Asteroid Institute is supported
through generous gifts from William K. Bowes Jr. Foundation, Tito’s
CHEERS, McGregor Girand Charitable Endowment, Galinsky Family,
D. Kaiser, P. Rawls Family Fund, J & M Montrym, Y. Wong and
Kimberly Algeri-Wong and three anonymous leadership donors plus
donors from 46 countries worldwide.
We would also like to thank Drs. Ed Lu, Aren Heinze, Mario Juric,´
and Siegried Eggl for taking the time and giving us invaluable feedback
during the expert interviews. Finally, we would also like to thank the
reviewers for this paper that have helped improve this manuscript.
This visualization system was implemented in the
OpenSpace system [6] and its source code is available at
https://github.com/OpenSpace/OpenSpace.
REFERENCES
[1] C. Acton, N. Bachman, B. Semenov, and E. Wright. A Look towards the
Future in the Handling of Space Science Mission Geometry. Planetary
and Space Science, 150:9–12, 2018. doi: 10.1016/j.pss.2017.02.013 5
[2] C. H. Acton. Ancillary Data Services of NASA’s Navigation and Ancillary
Information Facility. Planetary and Space Science, 44(1):65–70, 1996.
doi: 10.1016/0032-0633(95)00107-7 5
[3] A. Bischoff, M. Patzek, T. Di Rocco, A. Pack, A. Stojic, J. Berndt, and
S. Peters. Saint-Pierre-le-Viger (L5-6) from Asteroid 2023 CX1 Recovered
in the Normandy, France—220 years after the Historic Fall of L’Aigle
(L6 breccia) in the Neighborhood. Meteoritics & Planetary Science,
58(10):1385–1398, 2023. 8
[4] K. Bladin, E. Axelsson, E. Broberg, C. Emmart, P. Ljung, A. Bock, and
A. Ynnerman. Globe Browsing: Contextualized Spatio-Temporal Planetary Surface Visualization. IEEE Transactions on Visualization and
Computer Graphics, 21(1):802–811, 2018. doi: 10.1109/TVCG.2017.
2743958 9
[5] A. Boattini, G. D’abramo, G. Forti, and R. Gal. The Arcetri NEO precovery
program. Astronomy & Astrophysics, 375(1):293–307, 2001. 3
[6] A. Bock, E. Axelsson, J. Costa, G. Payne, M. Acinapura, V. Trakinski,
C. Emmart, C. Silva, C. Hansen, and A. Ynnerman. OpenSpace: A System
for Astrographics. IEEE Transactions on Visualization and Computer
Graphics, 2019. doi: 10.1109/TVCG.2019.2934259 2, 10
[7] A. Bock, A. Pembroke, M. L. Mays, L. Rastaetter, T. Ropinski, and
A. Ynnerman. Visual Verification of Space Weather Ensemble Simulations.
In 2015 IEEE Scientific Visualization Conference (SciVis), pp. 17–24, 2015.
doi: 10.1109/SciVis.2015.7429487 2
[8] N. Bowman and M. T. Heath. Computing Minimum-Volume Enclosing
Ellipsoids. Mathematical Programming Computation, 15:621–650, 2023.
doi: 10.1007/s12532-023-00242-8 5
[9] S. R. Chesley. Potential Impact Detection for Near-Earth Asteroids:
the Case of 99942 Apophis (2004 MN4). Proceedings of the International Astronomical Union, 1(S229):215–228, 2005. doi: 10.1017/
S1743921305006769 7
[10] V. Davydov. Tunguska Coals, Siberian Sills and the Permian-Triassic
Extinction. Earth-Science Reviews, 212:103438, 2021. 1
[11] F. DeMeo, C. Alexander, K. Walsh, C. Chapman, R. Binzel, et al. The
Compositional Structure of the Asteroid Belt. Asteroids iv, 1:13, 2015. 1
[12] F. Ferstl, K. Bürger, and R. Westermann. Streamline variability plots for
characterizing the uncertainty in vector field ensembles. IEEE Transactions on Visualization and Computer Graphics, 22(1):767–776, 2016. doi:
10.1109/TVCG.2015.2467204 2
[13] M. Granvik, J. Virtanen, D. Oszkiewicz, and K. Muinonen. OpenOrb:
Open-source Asteroid Orbit Computation Software including Statistical
Ranging. Meteoritics and Planetary Science, 44(12):1853–1861, 2009.
doi: 10.1111/j.1945-5100.2009.tb01994.x 4
[14] H. Guo, X. Yuan, J. Huang, and X. Zhu. Coupled Ensemble Flow Line
Advection and Analysis. IEEE Transactions on Visualization and Computer Graphics, 19(12):2733–2742, 2013. doi: 10.1109/TVCG.2013.144
2
[15] Ž. Ivezic, S. M. Kahn, J. A. Tyson, B. Abel, E. Acosta, et al. LSST: From ´
Science Drivers to Reference Design and Anticipated Data Products. The
Astrophysical Journal, 873(2):111, 2019. doi: 10.3847/1538-4357/ab042c
3
[16] R. L. Jones, C. T. Slater, J. Moeyens, L. Allen, T. Axelrod, K. Cook,
Ž. Ivezic, M. Juri ´ c, J. Myers, and C. E. Petry. The Large Synoptic Survey ´
Telescope as a Near-Earth Object Discovery Machine. Icarus, 303:181–
202, 2018. doi: 10.1016/j.icarus.2017.11.033 3
[17] A. Kamal, P. Dhakal, A. Y. Javaid, V. K. Devabhaktuni, D. Kaur, J. Zaientz,
and R. Marinier. Recent Advances and Challenges in Uncertainty Visualization: a Survey. Journal of Visualization, 2021. doi: 10.1007/s12650
-021-00755-1 2
[18] C. Kyba, S. Garz, H. Kuechly, A. S. De Miguel, J. Zamorano, J. Fischer,
and F. Hölker. High-resolution Imagery of Earth at Night: New Sources,
Opportunities and Challenges. Remote sensing, 7(1):1–23, 2015. 6
[19] F. Lan, M. Young, L. Anderson, A. Ynnerman, A. Bock, M. A. Borkin,
A. G. Forbes, J. A. Kollmeier, and B. Wang. Visualization in Astrophysics:
Developing New Methods, Discovering Our Universe, and Educating the
Earth. Computer Graphics Forum, 40(3):635–663, 2021. doi: 10.1111/cgf
.14332 2
[20] S. Larson, E. Beshore, R. Hill, E. Christensen, D. McLean, S. Kolar,
R. McNaught, and G. Garradd. The CSS and SSS NEO surveys. In
AAS/Division for Planetary Sciences Meeting Abstracts #35, vol. 35 of
AAS/Division for Planetary Sciences Meeting Abstracts, p. 36.04, 2003. 3
[21] S. Li, H. Cai, and V. R. Kamat. Uncertainty-aware Geospatial System for
Mapping and Visualizing Underground Utilities. Automation in Construction, 53:105–119, 2015. doi: 10.1016/j.autcon.2015.03.011 2
[22] L. Liu, L. Padilla, S. H. Creem-Regehr, and D. H. House. Visualizing
Uncertain Tropical Cyclone Predictions using Representative Samples
from Ensembles of Forecast Tracks. IEEE Transactions on Visualization
and Computer Graphics, 25(1):882–891, 2019. doi: 10.1109/TVCG.2018.
2865193 2
[23] M. Mirzargar, R. T. Whitaker, and R. M. Kirby. Curve Boxplot: Generalization of Boxplot for Ensembles of Curves. IEEE Transactions on
Visualization and Computer Graphics (TVCG), 20(12):2654–2663, 2014.
doi: 10.1109/TVCG.2014.2346455 2
[24] C. Norlund, H. Lewis, P. Atkinson, and J. Guo. NEOMiSS: A global
human vulnerability estimator and evacuation simulator. 3
[25] S. W. Paek, O. de Weck, J. Hoffman, R. Binzel, and D. Miller. Optimization and Decision-making Framework for Multi-Staged Asteroid
Deflection Campaigns Under Epistemic Uncertainties. Acta Astronautica,
167:23–41, 2020. doi: 10.1016/j.actaastro.2019.10.042 2
[26] S. W. Paek, O. de Weck, R. Polany, and P. Egger. Asteroid Deflection Campaign Design integrating Epistemic Uncertainties. In 2016
IEEE Aerospace Conference, pp. 1–14, 2016. doi: 10.1109/AERO.2016.
7500905 2
[27] A. Pang. Visualizing Uncertainty in Natural Hazards. In A. Bostrom,
S. French, and S. Gottlieb, eds., Risk Assessment, Modeling and Decision
Support: Strategic Directions, pp. 261–294. Springer Berlin Heidelberg,
Berlin, Heidelberg, 2008. doi: 10.1007/978-3-540-71158-2_12 2
[28] A. T. Pang, C. M. Wittenbrink, S. K. Lodha, et al. Approaches to Uncertainty Visualization. The Visual Computer, 13(8):370–390, 1997. doi: 10.
1007/s003710050111 2
[29] K. Potter, P. Rosen, and C. R. Johnson. From Quantification to Visualization: a Taxonomy of Uncertainty Visualization Approaches. In A. M.
Dienstfrey and R. F. Boisvert, eds., Uncertainty Quantification in Scientific
Computing, pp. 226–249. Springer Berlin Heidelberg, Berlin, Heidelberg,
2012. doi: 10.1007/978-3-642-32677-6_15 2
[30] C. M. Rumpf, H. G. Lewis, and P. M. Atkinson. Asteroid Impact Effects and Their Immediate Hazards for Human Populations. Geophysical
Research Letters, 44(8):3433–3440, 2017. doi: 10.1002/2017GL073191 1
[31] C. M. Rumpf, D. Mathias, D. Farnocchia, and S. Chesley. Projecting
Asteroid Impact Corridors onto the Earth. In 2019 IEEE Aerospace
Conference, pp. 1–7, 2019. doi: 10.1109/AERO.2019.8741912 3
[32] C. M. Rumpf, D. L. Mathias, L. F. Wheeler, J. L. Dotson, B. Barbee, J. Roa,
P. Chodas, and D. Farnocchia. Deflection Driven Evolution of Asteroid
Impact Risk Under Large Uncertainties. Acta Astronautica, 176:276–286,
2020. doi: 10.1016/j.actaastro.2020.05.026 3
[33] P. Schulte, L. Alegret, I. Arenillas, J. A. Arz, P. J. Barton, et al. The Chicxulub Asteroid Impact and Mass Extinction at the Cretaceous-Paleogene
Boundary. Science, 327(5970):1214–1218, 2010. doi: 10.1126/science.
1177265 1
[34] B. Shneiderman. The Eyes have it: A Task by Data Type Taxonomy for
Information Visualizations. In The craft of information visualization, pp.
364–371. Elsevier, 2003. 5
[35] D. Vokrouhlicky, A. Milani, and S. R. Chesley. Yarkovsky Effect on Small
Near-Earth Asteroids: Mathematical Formulation and Examples. Icarus,
148(1):118–138, Nov. 2000. doi: 10.1006/icar.2000.6469 3
[36] J. Wang, S. Hazarika, C. Li, and H.-W. Shen. Visualization and Visual
Analysis of Ensemble Data: a Survey. IEEE Transactions on Visualization
and Computer Graphics (TVCG), 25(9):2853–2872, 2019. doi: 10.1109/
TVCG.2018.2853721 2
[37] L. Wheeler, D. Mathias, J. L. Dotson, and NASA Asteroid Threat Assessment Project. Sensitivity of Asteroid Impact Risk to Uncertainty in
Asteroid Properties and Entry Parameters. In AAS/Division for Planetary
Sciences Meeting Abstracts #49, vol. 49 of AAS/Division for Planetary
Sciences Meeting Abstracts, p. 110.18, 2017. 2
[38] R. T. Whitaker, M. Mirzargar, and R. M. Kirby. Contour boxplots:
a method for characterizing uncertainty in feature sets from simulation ensembles. Transactions on Visualization and Computer Graphics,
19(12):2713–2722, 2013. doi: 10.1109/TVCG.2013.143 2
[39] D. Yeomans. Deflecting a Hazardous Near-Earth Object. Planetary
Defense Conference: Protecting Earth from Asteroids, 2009. 1, 2
[40] M. Zhang, L. Chen, Q. Li, X. Yuan, and J. Yong. Uncertainty-oriented Ensemble Data Visualization and Exploration Using Variable Spatial Spreading. IEEE Transactions on Visualization and Computer Graphics (TVCG),
27(2):1808–1818, 2021. doi: 10.1109/TVCG.2020.3030377 2
[41] S. Zhang, C. Demiralp, and D. H. Laidlaw. Visualizing diffusion tensor
mr images using streamtubes and streamsurfaces. IEEE Transactions
on Visualization and Computer Graphics, 9(4):454–462, 2003. doi: 10.
1109/TVCG.2003.1260740 2
</NEOviz: Uncertainty-Driven Visual Analysis of Asteroid Trajectories>


<ADVANCED MONTE CARLO SAMPLING TECHNIQUES FOR ORBITAL
CONJUNCTIONS ANALYSIS AND NEAR EARTH OBJECTS IMPACT PROBABILITY
COMPUTATION>
ADVANCED MONTE CARLO SAMPLING TECHNIQUES FOR ORBITAL
CONJUNCTIONS ANALYSIS AND NEAR EARTH OBJECTS IMPACT PROBABILITY
COMPUTATION
M. Losacco(1), M. Romano(1), P. Di Lizia(1), C. Colombo(1), R. Armellin(2), A. Morselli(3), and J.M. Sanchez P ´ erez ´
(3)
(1)Department of Aerospace Science and Technology, Politecnico di Milano, 20156, Milano, Italy, E-mail:
{matteo.losacco, matteo1.romano, pierluigi.dilizia, camilla.colombo}@polimi.it
(2)Surrey Space Centre, University of Surrey, GU2 7XH, Guildford, United Kingdom, E-mail: r.armellin@surrey.ac.uk
(3)European Space Operations Centre (ESOC/ESA), 64293, Darmstadt, Germany, E-mail: {alessandro.morselli,
jose.manuel.sanchez.perez}esa.int
ABSTRACT
This paper presents a survey of past and new results of
the application of advanced sampling techniques for orbital conjunction analysis and Near Earth Objects impact
probability computation. The theoretical background of
the methods is presented, along with the results of their
applications and a critical discussion of the benefits introduced.
Keywords: Space Debris; Conjunction Analysis; Near
Earth Objects; Impact Probability; Monte Carlo.
1. INTRODUCTION
When dealing with space debris and Near Earth Objects
(NEO), the major threat for human activities in space and
global safety is represented by possible collisions. The
risk of in-orbit collisions between operative satellites and
space debris, indeed, is a crucial issue in satellite operation. When a close approach is identified, it is necessary
to define an indicator that can reasonably tell how risky
the predicted conjunction may be, and a common practice for space agencies and satellite operators is to consider the collision probability for this purpose, together
with conjunction geometry and miss-distance [15]. On
the other hand, whenever a new asteroid is discovered, it
is of crucial importance to have accurate tools for detection and prediction of possible impacts. The task introduces relevant challenges due to the imperative of early
detection and propagation of its state and associated uncertainty, and the problem is made even more complicated by the fact that the dynamics describing the motion of these objects is highly nonlinear, especially during
close encounters with major bodies.
Present day approaches for orbital conjunctions analysis
and robust detection and prediction of planetary encounters and potential impacts by NEO include simplified
models or full nonlinear approaches. Simplified models significantly reduce the required computational burden, but their range of application is very narrow. On the
other hand, full nonlinear models, i.e. canonical Monte
Carlo (MC) simulations, represent a general and flexible way of approaching collision and impact probability
estimation, but they involve unavoidable intensive computations, which may be not suitable for satellite-debris
daily collision probability computation or fast NEO impact probability estimation. An elegant and effective
compromise was introduced by Milani, based on the concept of the Line Of Variations (LOV) [18].
In recent times, new advanced Monte Carlo techniques,
such as the Importance Sampling (IS), Line Sampling
(LS) and Subsets Simulation (SS) methods were developed, aiming at reducing the computational burden by
either restricting the sampling phase space or identifying optimal sampling paths within it. In parallel, an intense study on alternative uncertainty propagation methods was carried out, and innovative approaches based on
the use of differential algebra (DA) were proposed. This
paper presents a series of combinations and applications
of the mentioned techniques to orbital conjunctions analysis and NEO impact probability computation. The problem of orbital conjunctions is faced combining the mentioned SS and LS techniques with DA, obtaining a highorder, efficient tool that can be applied to any set of initial
states and using any arbitrary reference frame. The challenging task of NEO impact probability computation is
faced both directly applying the SS and LS techniques
and combining the IS method with DA. The presented
approaches offer competitive results in both fields, providing a significant improvement with respect to standard
MC performance while maintaining a sufficiently high
level of accuracy. For all the approaches, test cases are
presented, and a detailed analysis of the performance of
the proposed methods is offered.
The paper is organized as follows. First, a detailed theoretical description of differential algebra, Importance
Sampling, Line Sampling and Subset Simulation is ofProc. 1st NEO and Debris Detection Conference, Darmstadt, Germany, 22-24 January 2019, published by the ESA Space Safety Programme Office
Ed. T. Flohrer, R. Jehn, F. Schmitz (http://neo-sst-conference.sdo.esoc.esa.int, January 2019)
a, b ∈ R a, b ∈ F P
a × b
× ⊗
a ⊗ b
T
T
f, g
f × g
× ⊗
T
T
F, G
F ⊗ G
Figure 1: Analogy between the FP representation of real
numbers in computer environment (left) and the algebra
of Taylor polynomials in DA framework (right)[12].
fered. Then, the results of their application to orbital
conjunction analysis and NEO impact probability computation are presented, along with a critical discussion of
benefits and limitations of the presented approaches.
2. THEORETICAL BACKGROUND
A detailed description of the methods presented in the article is offered in this section. The differential algebra
approach and its recent enhancements are presented in
the first part. Then, three advanced orbital sampling techniques are described: Importance Sampling, Line Sampling and Subset Simulation.
2.1. Differential Algebra and Automatic Domain
Splitting
Differential algebra provides the tools to compute the
derivatives of functions within a computer environment
[25, 26, 23, 24, 16, 7]. The basic idea of DA is to bring the
treatment of functions and the operations on them to the
computer environment in a similar way as the treatment
of real numbers [7]. Real numbers, indeed, are approximated by floating point (FP) numbers with a finite number of digits. In a similar way, suppose two sufficiently
regular functions f and g are given. In the framework
of DA, these functions are converted into their Taylor series expansions, F and G respectively (see Fig. 1). For
each operation in the function space, an adjoint operation in the space of Taylor polynomials is defined. As a
results, the implementation of DA in a computer environment provides the Taylor coefficients of any function of
v variables up to a specific order n: any function f of v
variables can be expanded into its Taylor expansion up to
an arbitrary order n, along with the function evaluation,
with a limited amount of effort. The Taylor coefficients
of order n for sum and product of functions, as well as
scalar products with real numbers, can be directly computed from those of summands and factors. In addition to
basic algebraic operations, differentiation and integration
can be easily introduced in the algebra, thus finalizing
the definition of the differential algebra structure of DA
[5, 6]. The DA used in this work is implemented in the
DACE software [22].
A relevant application of DA is the automatic high or-
𝑡𝑡 = 𝑡𝑡0 𝑡𝑡 = 𝑡𝑡𝑖𝑖
error < 𝜀𝜀
𝑡𝑡 = 𝑡𝑡𝑖𝑖+1
error > 𝜀𝜀
error < 𝜀𝜀
Figure 2: ADS algorithm schematic illustration [28].
der expansion of the solution of an Ordinary Differential Equation (ODE) with respect to the initial conditions
[7, 13, 22]. This expansion can be achieved by considering that any integration scheme, explicit or implicit,
is characterized by a finite number of algebraic operations, involving the evaluation of the ODE right hand
side (RHS) at several integration points. By replacing
the operations between real numbers with those on DA
numbers, the nth order Taylor expansion of the flow of
the ODE, φ(t; δx0, t0) = Mφ(δx0), at each integration
time, for any perturbed initial condition x0 + δx0 can be
obtained.
While a single polynomial map can accurately describe
the evolution of an uncertainty set for short term propagation or dynamics characterized by a low level of nonlinearities, which is a common case of orbital conjunctions,
the accuracy of a single polynomial maps drastically decreases as the propagation window increases or the effect of nonlinearities becomes more important. This is
the typical case of NEO long term propagations. The approximation error of the DA representation is strictly related to the size of the domain the polynomial is defined
in [28]. This means that, if one divides the initial domain into smaller domains and then computes the Taylor
expansion around the center points of the new domains,
the error can be reduced, still covering the entire initial
set with all the generated polynomial expansions. Starting from these considerations, Automatic Domain Splitting (ADS) employs an automatic algorithm to determine
at which time ti
the flow expansion over the set of initial conditions is no longer able to describe the dynamics with enough accuracy [28]. Once this event has been
detected, the domain of the original polynomial expansion is divided along one of the expansion variables into
two domains of half their original size. By re-expanding
the polynomials around the new centre points, two separate polynomial expansions are obtained. In this way, the
splitting procedure guarantees a more accurate description of the whole uncertainty set at the current time epoch
ti
. After such a split occurs, the integration process is resumed on both generated subsets, until new splits are required. A representation of the ADS procedure is shown
in Fig. 2. The main degrees of freedom of the algorithm
are the tolerance for the splitting procedure and the maximum number of splits per subset. A detailed description
of the role of these parameters and ADS is offered in [28]
and [17].
2.2. Advanced orbital sampling techniques
Three advanced Monte Carlo sampling techniques are
presented in this section: Importance Sampling, Line
Sampling and Subset Simulation. All presented methods
aim at identifying the failure region of an uncertainty set
with the minimum amount of samples. For the problem
under study, failure is represented either by a collision in
orbital conjunction analysis or an impact for NEO impact
monitoring. The occurrence of the events is expressed by
the value of the performance function
gx(x0) = d(x0) − D (1)
where x0 represents a single initial condition for NEO
impact monitoring, a couple x0 = (x
1
0
, x
2
0
) for conjunction analysis, d is a function that maps an initial condition
x0 to a performance index, and D is the selected critical distance. For conjunction analysis, the performance
index is the Distance of Closest Approach (DCA). For
impact monitoring, it is the minimum planetocentric distance. The critical distance is a threshold beyond which
a “failure” is identified. For conjunction analysis, this
threshold can be identified as half of the diameter of the
sphere including both objects [20]. For impact monitoring, D represents the planet radius. According to this
definition, it follows that
gx(x0)



< 0 → failure
= 0 → limit state
> 0 → no failure
(2)
2.2.1. Importance Sampling
The Importance Sampling (IS) method amounts to replacing the original probability density function (pdf)
qx(x) with an Importance Sampling Distribution (ISD)
q˜x(x) arbitrarily chosen by the analyst so as to generate a large number of samples in the importance region
of the phase space F. This region represents the set of
initial conditions leading to collisions in orbital conjunction analysis or impacts for NEO impact monitoring. The
general IS algorithm is the following:
1. Identify a proper q˜x(x).
2. Express the failure probability P(F) as a function
of q˜x(x).
P(F) = Z
IF(x)qx(x)dx =
=
Z
IF(x)qx(x)
q˜x(x)
q˜x(x)dx
(3)
where IF(x) : IRv → {0, 1} is an indicator function
such that IF(x) = 1 if x ∈ F, 0 otherwise.
3. Draw NT samples x
k
: k = 1, 2, . . . , NT from the
importance sampling distribution q˜x(x). If a good
choice for the auxiliary pdf is made, the generated
samples concentrate in the region F.
4. Compute the estimate Pˆ(F) for the failure probability P(F) by resorting to equation (3):
Pˆ(F) = 1
NT
X
NT
k=1
IF(x
k
)qx(x
k
)
q˜x(xk)
(4)
5. Compute the variance of the estimator Pˆ(F) as:
σ
2
(Pˆ(F)) = 1
NT
Z
I
2
F
(x)q
2
x(x)
q˜
2
x(x)
q˜x(x)dx
− P
2
(F)

≈
1
NT

Pd2(F) − Pˆ2
(F)

(5)
The selection of the ISD represents the most critical point
for the method. Several techniques have been developed
in order to find the one giving small variance for the estimator [30]. As later described in Sec. 4, this sampling
method can be combined with DA for NEO impact monitoring, and the ISD can be modelled according to the
results of the DA propagation.
2.2.2. Line Sampling
The Line Sampling (LS) method is a Monte Carlo-based
approach for the estimation of small probabilities. The
main idea behind LS is to transform a high-dimensional
problem into a number of conditional one-dimensional
problems solved along an important direction α. The LS
method follows four steps: 1) the mapping of random
samples from the physical coordinate space into a normalized standard space, 2) the determination of the reference direction α, 3) the probing of the failure region
along lines following the reference direction α, 4) the estimation of the failure probability.
In the LS approach, a vector of uncertain parameters
x0 ∈ IRn
is first transformed into the adjoint vector
θ ∈ IRn
belonging to the so-called “standard normal
space”, where each variable is represented by an independent central unit Gaussian distribution. This is done
by relying on Rosenblatt’s transformation [27]
θ = Tx,θ(x0)
x0 = Tθ,x(θ)
(6)
A natural choice for the important direction is the normalized gradient of the performance function at the nominal
point in the standard normal space
α =
∇θgθ(θ)
||∇θgθ(θ)||2
(7)
If not available analytically, the gradient can be estimated
numerically. Another approach consists in obtaining an
estimate by computing the normalized centre of mass of
the failure domain. This is achieved by Monte Carlo
Markov Chain (MCMC), using as seed a point belonging to the failure region or close to it, and computing the
mean of the Ns samples generated in the failure region.
Once the important direction is identified, the LS method
proceeds as follows [20]:
• Sample NT vectors θ from the normal multidimensional joint probability distribution.
• For each sample i, estimate its conditional onedimensional failure probability Pˆi
. This task is
solved performing the following operations (see also
Fig. 3)
– Project the vector θ
i onto the straight line passing through the origin and perpendicular to α
to obtain the vector θ
i,⊥.
– Write the parametric equation of samples
along the important direction θ
i = θ
i,⊥ + c α
– Compute the values of c
i
j
, with j = 1, 2,
for which the performance function is equal
to zero. This step requires the evaluation of
the performance function, which involves extra numerical propagations.
– If the two values are coincident or no solution
is found, then the ith one-dimensional probability Pˆi
is equal to zero, otherwise, given
the two solutions c
i
1
and c
i
2
, the probability becomes
Pˆi
(F) = P[c
i
2 ≤ N(0, 1) ≤ c
i
1
] =
= Φ(c
i
1
) − Φ(c
i
2
)
where Φ(c
i
j
) is the standard normal cumulative
distribution function, N(0, 1) is the standard
normal distribution, with zero mean and unit
standard deviation, and F is the failure event
(in-orbit collision for conjunction analysis or
NEO impact for impact monitoring).
• Compute the unbiased estimator PˆNT (F) as the
sample average of the independent conditional onedimensional probability estimate
PˆNT
(F) = 1
NT
X
NT
i=1
Pˆi
(F) (8)
The variance of the estimator is given by
σ
2
(PˆNT
(F)) = 1
NT(NT − 1)
X
NT
i=1
(Pˆi
(F)−PˆNT
(F))2
(9)
2.2.3. Subset Simulation
The Subset Simulation (SS) method is a Monte Carlo
method based on the principle of computing small failure
Figure 3: Scheme of the iterative procedure used to sample each line in the standard normal coordinate space.
The failure region is labeled with F, with a single border
highlighted as a red line (image courtesy of [31]).
probabilities as the product of larger conditional probabilities [3, 8, 32]. Given a target failure event F, let
F1 ⊃ F2 ⊃ ... ⊃ Fn = F be a sequence of intermediate failure events, so that Fk = ∩
k
i=1Fi
. Considering
a sequence of conditional probabilities, then the failure
probability can be written as
P(Fn) = P(F1)
nY−1
i=1
P(Fi+1|Fi) (10)
where P(Fi+1|Fi) represents the probability of Fi+1
conditional to Fi
.
The method is initialized using standard Monte Carlo
to generate samples at the so called conditional level 0
(CL 0) starting from the available nominal state vectors
and related uncertainties of the investigated objects. The
number of samples generated at this level is maintained
for each generated conditional level and it is referred to
as N. Once the failure region F1 is identified, a MCMC
Metropolis Hastings algorithm is used to generate conditional samples in the identified intermediate failure region. Another intermediate region is then located, and
other samples are generated by means of MCMC. The
procedure is repeated until the failure region is identified. A scheme of the method is shown in Fig 4. In the
presented approach, the intermediate failure regions are
identified by assuming a fixed value of conditional probability P(Fi+1|Fi) = p0. The identification of each conditional level, therefore, is strictly related to this value,
and changes accordingly step by step, as explained in the
followings. The resulting Subset Simulation algorithm
follows the general description offered in [20] and goes
through the following steps:
1. Set i = 0 and generate N samples x
0,k
0
, k =
0, ..., N at conditional level 0 by standard Monte
Carlo.
a)
Matteo Losacco, PhD annual evaluation 1
𝑥𝑥1
𝑥𝑥2
Standard MC
Impact region 𝐹𝐹
b)
Matteo Losacco, PhD annual evaluation 1
𝑥𝑥1
𝑥𝑥2
𝑃𝑃(𝐹𝐹1) = 𝑝𝑝0
𝐹𝐹1
Impact region 𝐹𝐹
Standard MC
c)
Matteo Losacco, PhD annual evaluation 1
𝑥𝑥1
𝑥𝑥2
𝑃𝑃(𝐹𝐹1) = 𝑝𝑝0
Impact region 𝐹𝐹
𝐹𝐹1
Standard MC
d)
Matteo Losacco, PhD annual evaluation 1
𝑥𝑥1
𝑥𝑥2
𝑃𝑃(𝐹𝐹1) = 𝑝𝑝0
Impact region 𝐹𝐹
𝐹𝐹1
𝐹𝐹2
𝑃𝑃(𝐹𝐹2) = 𝑝𝑝0
𝑃𝑃(𝐹𝐹3) = 𝑁𝑁𝑛𝑛/𝑁𝑁
Standard MC
Figure 4: Subset Simulation process: a), initialization by standard MC, b), CL 1 identification, c), samples generation by
means of MCMC, d), new iterations and impact region identification (image courtesy of [17]).
2. Propagate each sample for the selected time window
and compute its performance function, defined as
g
i
x(x0) = d(x0) − Di

< 0 → x0 ∈ (i)th CL
> 0 → x0 ∈/ (i)th CL
(11)
where Di
is the current threshold, whereas d is a
function that provides, for each sample, its performance index. For conjunction analysis, the index
is the Distance of Closest Approach, for NEO impact probability computation, this index represents
the minimum planetocentric distance.
3. Sort the N samples in descending order according
to their associated value of performance function.
4. Identify an intermediate threshold value Di+1 as the
(1−p0)Nth element of the list. Define the (i+ 1)th
conditional level as Fi+1 = {d < Di+1}. Considering how the threshold was defined, the associated
conditional probability P(Fi+1|Fi) = p0
5. If Di+1 < D, go the last step, otherwise select the
last p0N samples of the list x
i,j
0
, j = 0, ..., p0N.
By definition, these samples belong to the (i + 1)th
conditional level.
6. Using MCMC, generate (1−p0)N additional conditional samples starting from the previously selected
seeds belonging to Fi+1. A sample is set to belong
to Fi+1 according to its performance function value.
7. Set i = i + 1 and return to step 2
8. Stop the algorithm.
The total number of generated samples is
NT = N + (n − 1)(1 − p0)N (12)
where n is the overall number of conditional levels required to reach the failure region. Since the conditional
probability is equal to p0 for each level, the failure probability becomes:
P(F) = P(Fn) = P(Fn|Fn−1)p
n−1
0 =
= p
n−1
0 NF/N
(13)
where NF is the number of samples belonging to the failure region.
A Bayesian post-processor for SS (SS+) is suggested in
[32] to refine the computed failure probability and determine higher moments. If we define
nl =

p0N if l < n
NF if l = n
(14)
Table 1: Time, distance and velocity of closest approach for the Keplerian test cases and reference values for collision
probability [20].
Test case TCA (days) DCA (m) ∆vTCA (m/s) D (m) Pˆ
c (MC) NT [11] Pˆ
c [2] % err (-)
5 2.0 2.449 0.520 10 4.454e-2 2.30e6 4.440e-2 -0.32%
6 2.0 2.449 0.173 10 4.340e-3 2.50e7 4.324e-3 -0.36%
7 2.0 3.183 0.196 10 1.614e-4 6.71e8 1.580e-4 -2.13%
the first moment of the distribution of the failure probability becomes
ESS+ [P] = Yn
l=1
nl + 1
N + 2
(15)
whereas the second moment is expressed by
ESS+ [P
2
] = Yn
l=1
(nl + 1) (nl + 2)
(N + 2)(N + 3) (16)
The variance of the estimator, therefore, can be computed
as
σ
2
(P) = E[P
2
] − (E[P])2
(17)
Equations 15 and 17 will represent our reference for the
analyses presented in this paper.
3. ORBITAL CONJUNCTION ANALYSIS
The estimation of the collision probability for in-orbit
conjunction analysis requires the computation of a multivariate integral. Different methods exist for the computation of this multi-dimensional integral[1, 4, 21]. Most of
these approaches rely on strong assumptions, assuming
uncorrelated position uncertainties for the two objects,
objects moving along straight lines at constant velocity
during the conjunction, negligible uncertainties in velocity, and constant, multi-dimensional Gaussian uncertainties in position during the whole encounter for both objects. These assumptions produce accurate results when
the relative motion between the satellite and the object is
rectilinear and the conjunction occurs close to the initial
epoch.
Methods that account for nonlinearities, which are usually relevant for GEO conjunctions, typically rely on analytical methods or Monte Carlo simulations. Despite being a general and flexible way to compute collision probability, the MC approach has the main drawback of requiring intensive computation. For this reason, MC methods
are not suitable for daily collision probability computation, since the results can be obtained in a timely manner
only relying on simplified orbital dynamics, such as two
body propagators or SGP4/SDP4.
The approach illustrated in this paper was first presented in [20] and combines the previously described
DA method with standard MC, LS and SS, obtaining the
DAMC, DALS and DASS methods for orbital conjunction analysis. The basic idea is to use DA to obtain the
Taylor expansion of the Time of Closest Approach (TCA)
and the Distance of Closest Approach (DCA) of the two
orbiting objects. The occurrence of close approaches
is first identified using the nominal initial orbital states.
Then, DA is used to propagate sets of initial conditions
by computing the Taylor expansion of the final states x
1
f
and x
2
f
at the nominal TCA. The resulting polynomials
are functions of both the final time and the initial uncertain state vectors x
1
0
and x
2
0
. The polynomial map of
the relative distance between the two objects is then computed through simple algebraic manipulations. By using
partial polynomial inversion techniques and imposing the
stationarity condition of the relative distance with respect
to time, the Taylor expansions of TCA and DCA with respect to x
1
0
and x
2
0
.
[t
∗
] = t
∗ + Mt
∗ (δx
1
0
, δx
2
0
)
[d
∗
] = d
∗ + Md∗ (δx
1
0
, δx
2
0
)
(18)
are computed. Given any perturbed initial condition of
the two objects (δx
1
0
, δx
2
0
), the evaluation of the two
Taylor polynomials provides the values of the TCA and
DCA. This enables a drastic reduction of the computational cost by replacing standard propagations with multiple polynomial evaluations. At this stage, the sampling is
performed relying either on MC, LS or SS, taking advantage of the availability of the resulting polynomial maps.
In particular, for the DALS method, the availability of
the DCA expansion is exploited for the computation of
the important direction α, whereas all numerical propagations of the standard algorithm are replaced by numerical evaluations. A detailed explanation is offered in [20].
For DAMC and DASS, the advantages are limited to the
replacement of numerical propagations.
We present here the results of the application of the
proposed methods for three different test cases taken
from [2], where simple Keplerian dynamics is used for
performance comparison. Results and figures are taken
from [20]. The comparison is done considering the number of propagated samples NT, the standard deviation of
the estimated collision probability σ, the coefficient of
variations δ = σ/Pˆ, and the figure of merit ∆, defined as
∆ = δ
p
NT (19)
The selected figure of merit does not depend on the number of samples, since for Monte Carlo methods σ ∝
Table 2: Computed collision probability for the Keplerian test cases. The comparison is shown in terms of relative error
with respect to reference Pˆ
c, number of samples NT, computational time tCPU, coefficient of variations δ and figure of
merit ∆ [20].
Test case Method Pˆ
c (-) % err (-) NT tCPU (s) δ (-) ∆ (-)
5
MC 4.452e-2 -0.05% 1.0e5 4.75 1.465e-2 4.63
DAMC 4.459e-2 +0.11% 1.0e5 0.67 1.464e-2 4.63
DALS 4.451e-2 -0.07% 5.0e3 2.53 7.662e-4 0.05
DASS 4.450e-2 -0.09% 2.0e4 0.13 2.738e-2 3.87
6
MC 4.339e-3 -0.01% 1.0e6 43.21 1.515e-2 15.14
DAMC 4.350e-3 +0.24% 1.0e6 6.67 1.513e-2 15.13
DALS 4.341e-3 +0.03% 5.0e3 2.58 1.484e-3 0.11
DASS 4.328e-3 -0.27% 4.0e4 0.27 3.586e-2 7.17
7
MC 1.615e-4 +0.04% 2.7e7 1155.36 1.514e-2 78.68
DAMC 1.612e-4 -0.15% 2.7e7 179.34 1.516e-2 78.76
DALS 1.621e-4 +0.41% 5.0e3 1.43 1.936e-2 1.37
DASS 1.626e-4 +0.72% 6.0e4 0.43 4.580e-2 11.22
particular, it is worth noting that for order k ¼ 4 it also
exceeds the computational time of pointwise MC, that
for test case 5 requires a lower number of samples. Nevertheless, the computational effort of DALS and DASS
decreases for lower collision probability, becoming nearly
103 times lower than the one of a standard Monte Carlo
method for test case 7.
The figure of merit D, normalized for each test case with
the value of the standard Monte Carlo method, is plotted
against the collision probability Pc in Fig. 8(a). The same
criteria used in Fig. 7 for markers coloring and shape is
used, i.e. different markers are used for each method and
they are colored according to the expansion order using a
gray scale. The normalized unitary c.o.v. is equal to 1 for
the DAMC since the same number of samples of the standard Monte Carlo method is used. The lower value is
achieved with DALS, which is two order of magnitude
lower than DAMC. The efficiency of DASS increases for
lower probabilities. The use of different expansion orders
does not affect the final value of the normalized D, since
points are overlapping and indistinguishable. The only
exception is found for DALS in cases 5 and 6, where the
normalized D is slightly higher when k ¼ 1. This is probably due to a slightly lower accuracy of the first-order
DCA expansion in those cases.
The figure of merit X is plotted versus the collision probability in Fig. 8(b), again normalized with respect to the
value obtained with the standard Monte Carlo method.
Since the computational time increases with the expansion
order as shown in Fig. 7, the value of X decreases for
Fig. 6. Comparison of the collision probability obtained with the tested methods for the Keplerian test cases. The DA expansion order is k ¼ 3 for
DAMC, DALS, and DASS. The solid black line is the reference value for collision probability and the gray lines are the 5% relative error bounds.
Fig. 7. Normalized computational time of DAMC, DALS, and DASS for
the Keplerian test cases vs. collision probability for different expansion
orders. Markers are colored using a grayscale, where black is used for the
expansion order k ¼ 1 and lighter gray for k ¼ 4.
Case 7 Case 6
Case 5
(a)
higher order. The best performing among the three
methods for the considered test cases is DALS, since the
normalized figure of merit is at least 10 times larger than
the one of DASS and 102 times larger than DAMC. For
test case 7, where collision probability is lower, the efficiency of DASS in terms of X is higher than the other cases.
Note that the value of X for an expansion order k ¼ 1 is
lower than the one obtained with k ¼ 2 for DALS for test
case 5 and 6. As stated before, the reason is the slightly
lower accuracy of the map in this case, that is not mitigated
by the lower computational time.
To conclude this analysis, the collision probability computed with Alfano’s formula and the DA-based methods
with an expansion order k ¼ 1 are compared in Table 3.
It can be observed that using a first-order DA expansion
the percentage relative error is similar to that obtained
using Alfano’s method for test case 7, where the relative
motion is no more linear. The DA-methods at first
order are therefore equivalent to Alfano’s analytical
approximation.
4.2. Comparison of the methods on real conjunctions
In this section four test cases are considered to test the
algorithms for collision probability computation. The
selected test caswith different rcases are listedconjunction calisted in the secreport the TCAthe collision thrOn the last cousing Alfano’sused for orbit pdix A.
The collisionDALS, and DAlisted in Table 5ity are considerinitial positionsobservations geare given in Apshown that, foof a third orderthe collision throne week.
The numberusing Eq. (10) cthe number of sFig. 8. Normalized figures of merit D and X of DAMC, DALS, and DASS for the Keplerian tesorders. Markers are colored using a gray scale, where black is used for the expansion order k ¼ 1Table3vances in Space Research 55 (2015) 311–333
(b)
Figure 5: omputational time tCP U (a) and figure of merit ∆ (b) of DAMC, DALS and DASS for the three
considered test cases vs collision probability for different expansion orders. Markers are colored using a grey scale, where
black is used for the expansion order k = 1 and light gray for k = 4 [20].
1/
√
NT. It is designed to enable the comparison of different methods in terms of accuracy and number of samples
required to reach that accuracy level. The three considered test cases are labeled as test case 5, 6 and 7, and are
characterized by linear relative motion between the two
objects, motion at the boundary of linear relative motion
and nonlinear relative motion respectively.
Table 1 shows the characteristics of the three selected
cases in terms of time, distance, and relative velocity
at the closest approach ∆vTCA, along with the reference
value for the collision probability as obtained with standard Monte Carlo simulation. The value of NT used for
MC simulations was selected according to [11]. For each
trial, two sets of initial conditions are sampled from each
initial covariance matrix and the associated DCA is computed in the proximity of the nominal TCA. The last two
columns of the table show the collision probability obtained using Alfano’s formula and its associated percentage relative error to the reference.
Table 2 shows the results obtained with the three proposed DA-based approaches. The results are expressed
in terms of expected collision probability Pˆ
c, relative error with respect to reference, number of required samples
NT, computational time tCPU, coefficient of variation δ
and figure of merit ∆. For all simulations, an expansion
order k for DA propagation equal to 3 was considered.
As can be seen, for all three test cases, the computed collision probability values are in good agreement with the
Table 3: Nominal equinoctial parameters and related uncertainty covariance matrix for asteroid 2017 RH16 at 6475
MJD2000 (September 24, 2017).
a (AU) P1 (-) P2 (-) Q1 (-) Q2 (-) l (deg)
State 0.8752 -0.1867 -0.4014 -0.0020 0.0050 319.9653
Cov.
9.7910e-08 1.2616e-08 3.6644e-08 -1.1067e-09 -2.2278e-10 -2.2745e-06
1.2616e-08 1.8745e-09 5.2329e-09 -1.3999e-10 -3.7340e-11 -2.1989e-07
3.6644e-08 5.2329e-09 1.4767e-08 -4.0882e-10 -1.0114e-10 -7.0073e-07
-1.1067e-09 -1.3999e-10 -4.0882e-10 1.2541e-11 2.4080e-12 2.6480e-08
-2.2278e-10 -3.7340e-11 -1.0114e-10 2.4080e-12 9.2142e-13 2.6346e-09
-2.2745e-06 -2.1989e-07 -7.0073e-07 2.6480e-08 2.6346e-09 7.4365e-05
reference values. A first consideration on the use of differential algebra can be made by comparing the performance of MC and DAMC for all three selected test cases.
By analyzing the required computational time, indeed,
it is evident that, as the expected probability decreases,
the required number of samples increases, and the advantages guaranteed by the introduction of the polynomial
map become more and more evident. A second consideration can be done by comparing the results of DAMC
with DALS and DASS. As can be seen, for all three test
cases, both methods guarantee good accuracy levels with
a number of samples that is at least one order of magnitude lower than what required by MC and DAMC. This
difference increases as the expected collision probability
decreases. This fact has obviously a direct positive consequence on the value of the figure of merit ∆.
The performance of all three DA-based methods is
strongly influenced by the selected expansion order. As
the expansion order increases, the accuracy obtained in
the definition of the polynomial maps increases, but the
required computational time increases as well. Figure 5a
shows the computational time required by the three DAbased methods, normalized with respect to the computational time of standard MC, as a function of the collision
probability, for expansion orders from 1 to 4. A greyscale
is used, with black for k = 1 and light gray for k = 4.
Different markers are used for the three methods: squares
for DAMC, circles for DALS and triangles for DASS.
Let us concentrate the analysis on the performance of
DALS and DASS. As can be seen, for each method, the
required computational time increases as the expansion
order increases, and the advantages with respect to standard MC increase as the expected collision probability
decreases. In particular, DALS becomes progressively
more and more appealing with decreasing collision probabilities.
Figure 5b shows the figure of merit ∆, normalized for
each test with the value of the standard Monte Carlo
method, as a function of the collision probability. The
same criteria used in Fig. 5a are used here. The normalized unitary coefficient of variations is equal to 1 for the
DAMC since the same number of samples of standard
MC is used. The lower value is achieved with DALS,
which is two orders of magnitude lower than DAMC. The
efficiency of DASS increases for lower probabilities. As
can be seen, the effect of the expansion order in this case
is negligible. Overall, both DALS and DASS offer a valuable alternative tool with respect to standard MC for conjunction analysis.
4. NEO IMPACT PROBABILITY COMPUTATION
Present day approaches for robust detection and prediction of planetary encounters and potential impacts by
NEO mainly refer to linearized models [9] or full nonlinear orbital sampling [29, 10, 19]. Among nonlinear methods, the Line of Variations (LOV) method [19] represents
the reference technique for impact probability computation of NEO.
The preferred approach to detecting potential impacts depends on the uncertainty in the estimated orbit, the investigated time window and the dynamics between the
observation epoch and the epoch of the expected impact
[14]. Linear methods are preferred when linear approximations are reliable for both the orbit determination and
uncertainty propagation. When these assumptions are not
valid, one must resort to more computationally intensive
techniques: among these, Monte Carlo methods are the
most accurate but also the most computationally intensive, whereas the LOV method guarantees compute times
3-4 orders of magnitude lower than those required in MC
simulations, though the LOV analysis may grow quite
complex after it has been stretched and folded by multiple close planetary encounters, leaving open the possibility of missing some pathological cases [14].
Two different approaches are described in this paper: for
the general problem of NEO impact probability computation, the presented LS and SS methods are directly applied, providing an efficient tool for impact monitoring.
The special case of resonant return analysis is faced using a dedicated combination of ADS and IS. For each
method, a single test case is presented. The propagations are carried out in Cartesian coordinates with respect to the J2000 reference frame centered in the Solar System Barycenter (SSB), with the inclusion of the
a)
No impact
Impact
b) c)
Figure 6: Visualization of the initial dispersion in the uncertainty space (δα,δl) in the case of asteroid 2017 RH16: a)
initial conditions leading to impact obtained via standard MC, b) boundaries of the subdomain identified via LS, c)
samples per conditional level obtained with SS.
Table 4: Application of Line Sampling and Subset Simulation to the case of asteroid 2017 RH16 against standard Monte
Carlo simulation.
NIC NT Pˆ (-) σˆ (-) δ (-) ∆ (-)
MC 5e4 5e4 1.42e-3 1.68e-4 1.19e-1 26.61
LS (ref) 1e3 8245 1.56e-3 7.19e-5 4.60e-2 4.18
LS (σ
MC) 164 1557 1.43e-3 1.70e-4 1.19e-1 4.70
LS (N MC
T
) 6250 5.02e4 1.52e-3 2.85e-5 1.89e-2 4.19
SS (ref) 1e3 2.8e3 1.45e-3 2.24e-4 1.55e-1 8.19
SS (σ
MC) 2e3 5.6e3 1.43e-3 1.57e-4 1.10e-1 8.27
SS (N MC
T
) 1.8e4 5.04e4 1.42e-3 5.20e-5 3.66e-2 8.22
gravitational contributions of the Sun, all the major planets, and the Moon, including relativity effects [28]. All
the physical constants (gravitational parameters, planetary radii, etc.) and ephemerides are obtained from the
JPL Horizons database via the SPICE toolkit (https:
//naif.jpl.nasa.gov/naif/). All propagations
are carried out in dimensionless units (with the scaling
length and time respectively equal to 1 AU and 1 solar
year) using the adaptive Dormand-Prince Runge-Kutta
scheme of 8th order (RK78), with absolute and relative
tolerances both equal to 1e-12.
The first test case presented is asteroid 2017 RH16. Asteroid 2017 RH16 will encounter the Earth at short distance in 2026, with a currently expected impact probability of 1:689. Table 3 shows the initial conditions in terms
−10 −8 −6 −4 −2 0 2 4 6 8 10
−4
−2
0
2
4
δa (AU)
δ
l (rad)
Date
(TDB)
2036−
05−31.0
2036−
03−17.0
2036−
04−01.0
2036−
05−01.0
2036−
05−16.0
2036−
03−02.0
2036−
04−16.0
x 10−8
x 10−6
Figure 7: Projection of the generated subsets onto the a − l plane of the initial conditions for test case of asteroid 99942
Apophis (ADP propagation, order 5, tolerance 1e-10, Nmax 10, 5σ domain). In light blue, boundaries of the ISD [17].
Table 5: Apophis equinoctial parameters and related uncertainties on June 18, 2009 00:00:00 (TDB).
Nominal value σ
a 0.922438242375914 2.29775e-8 AU
P1 -0.093144699837425 3.26033e-8 -
P2 0.166982492089134 7.05132e-8 -
Q1 -0.012032857685451 5.39528e-8 -
Q2 -0.026474053361345 1.83533e-8 -
l 88.3150906433494 6.39035e-5 ◦
of nominal equinoctial parameters and related uncertainties, as obtained from the Near Earth Objects Dynamic
Site1
. This case was used to test the presented LS and SS
method against the computation of the impact probability
value in 2026.
Figure 6 shows the results of the numerical simulations
obtained with standard MC, LS and SS, in terms of samples distribution on the semi-major axis- true longitude
(a − l) plane of initial conditions. Figure 6(a) shows the
results of standard MC, with impacting samples represented in red. Fig. 6(b) shows the results obtained with
LS. Grey dots represent samples drawn from the initial
distribution that do not lead to impact, whereas green
dots are the initial conditions along the boundaries of the
impact regions resulting from the LS application. Figure 6(c), finally, shows the evolution of the conditional
samples obtained with SS. The method was applied by
employing 1000 samples per conditional level, a fixed
conditional probability equal to 0.2, and an auxiliary distribution centered in the current sample and with the same
magnitude of the original one. Blue dots represent samples drawn at CL0 by standard MC, whereas different
1http://newton.dm.unipi.it/neodys/
colors are used for each conditional level. The algorithm
rapidly identifies the impact region, and after three conditional levels, the current threshold distance becomes
lower than the Earth radius, and the algorithm stops.
An overview of the obtained performance is shown in Table 4. For both LS and SS, three results are presented:
the results at convergence, the results obtained using a
number of samples granting the same accuracy level of
standard MC (σ
MC), and the results obtained performing
the same number of propagations of standard MC (N MC
T
).
The comparison between the methods is performed by
analyzing the number of random initial conditions NIC,
the total number of orbital propagations NT, the impact
probability estimate Pˆ, the standard deviation σˆ of Pˆ, the
coefficient of variation δ of the probability estimate and
the figure of merit ∆. The overall number of propagations NT was selected as a term of comparison for the
computational burden of the methods. As can be seen,
both methods allows us to improve the performance of
standard MC both in terms of achievable accuracy and
computational cost. That is, the same accuracy level of
standard MC simulation can be obtained by performing
a number of propagations that is one order of magnitude
lower for both LS and SS, or a higher accuracy level can
be obtained by propagating the same number of samples
used for standard MC. This result is confirmed by the
value of the figure of merit ∆.
The critical case of impact probability computation of
NEO during resonant returns is here faced combining
ADS and IS. The method, referred to as ADP-IS, was
presented in [17] and can be divided in two phases. During the first phase, the epoch of the resonant return is estimated with DA, and the uncertainty set is propagated
in time by means of ADS. A tailored pruning action is
introduced during the propagation, and only subsets that
are involved in the resonant return of interest are maintained throughout the propagation. At the end of the first
Table 6: Comparison between standard MC and ADP-IS method for the test case of the resonant return of asteroid 99942
Apophis [17].
NT tCPU / sample (s) Pˆ (-) σˆ (-)
MC 1.0e6 1.2 2.20e-5 4.71e-6
ADP-IS 341804 0.012 1.90e-5 6.81e-6
phase, several subdomains are generated only in specific
regions of the initial uncertainty set. At this point, the
sampling phase starts. An ISD is shaped on the basis of
the generated subsets, and the IS phase is started. Samples are propagated to the epoch of the expected resonant
return, and impacting samples are recorded. The procedure terminates when the variations in the estimated impact probability becomes lower than a selected threshold.
The method is here presented for the critical case of asteroid 99942 Apophis. Table 5 shows the nominal initial state and associated uncertainties σ for Apophis on
June 18, 2009 expressed in terms of equinoctial parameters p = (a, P1, P2, Q1, Q2, l), considering a diagonal
covariance matrix. Data were obtained from the Near
Earth Objects Dynamic Site in September 2009. Asteroid
Apophis will have a close encounter with Earth on April
13, 2029 with a nominal distance of 3.8e4 km. According to the selected initial conditions, though an impact in
2029 can be ruled out, the perturbations induced by the
encounter open the door to resonant returns in 2036 and
2037. The aim is therefore to apply the presented method
to provide an estimate for the impact probability at the
epoch of the first resonant return, in 2036.
Figure 7 shows the results of the first propagation phase
in terms of subdomains distribution on the a − l plane.
The method was applied considering an uncertainty set
of 5σ size, an expansion order of the ADS propagation
equal to 5, a maximum number of splits equal to 10, and a
tolerance for the splitting procedure equal to 1e-10. Each
subdomain is coloured according to the epoch at which
its DA propagation was stopped. The subdomains distribution represents the starting point for the second phase,
based on IS. First, an ISD including all the generated subsets is defined. Figure 7 shows in light blue the boundaries of this ISD as projected onto the a − l plane. Then,
the sampling phase starts: samples are generated on the
basis of the ISD, each sample is associated, if possible,
to a specific subset, and it is propagated to the epoch of
minimum geocentric distance. All impacting samples are
stored, and an estimate for the impact probability is obtained. Table 6 shows the results for the considered case,
and a comparison with standard MC. Given the different
propagation window for samples of the two method, an
equivalent computational time per sample is introduced
here as a term of comparison. As can be seen, the combination of ADS and IS guarantees a significant reduction
in the required computational time, still providing good
accuracy levels.
5. CONCLUSIONS
This paper presented an overview of past and present results obtained applying advanced uncertainty propagation
and sampling techniques to the challenging tasks of inorbit conjunction analysis and NEO impact probability
computation. All presented approaches offer great benefits in terms of computational effort with respect to standard Monte Carlo approach, while granting a competitive
accuracy level. Our research activity is currently devoted
to the enhancement of the presented NEO impact probability computation tools, with the aim of widening the
range of test cases and obtaining a more robust tool for
impact monitoring.
ACKNOWLEDGMENTS
M. Losacco and M. Romano gratefully acknowledge professors E. Zio, N. Pedroni and F. Cadini from Politecnico di Milano for their introduction to the Line Sampling
and Subset Simulation methods. The part of work about
the Line Sampling was funded by the European Research
Council (ERC) under the European Unions Horizon 2020
research and innovation programme (grant agreement No
679086 COMPASS).
REFERENCES
1. Akella, M. R. and Alfriend, K. T. (2000). Probability
of collision between space objects. Journal of Guidance, Control, and Dynamics, 23(5):769–772.
2. Alfano, S. (2009). Satellite conjunction Monte
Carlo analysis. Advances in Astronautical Sciences,
134:2007–2024.
3. Au, S.-K. and Beck, J. L. (2001). Estimation of
small failure probabilities in high dimensions by Subset Simulation. Probabilistic Engineering Mechanics,
16(4):263–277.
4. Berend, N. (1999). Estimation of the probability of
collision between two catalogued orbiting objects. Advances in Space Research, 23(1):243–247.
5. Berz, M. (1986). The new method of TPSD algebra for the description of beam dynamics to high orders. techreport AT-6:ATN-86-16, Los Alamos National Laboratory.
6. Berz, M. (1987). The method of power series tracking for the mathematical description of beam dynamics.
Nuclear Instruments & Methods in Physics Research,
Section A: Accelerators, Spectrometers, Detectors, and
Associated Equipment, A258(3):431–436.
7. Berz, M. (1999). Modern Map Methods in Particle Beam Physics. Advances in Imaging and Electron
Physics. Academic Press, 1 edition.
8. Cadini, F., Avram, D., Pedroni, N., and Zio, E. (2012).
Subset Simulation of a reliability model for radioactive
waste repository performance assessment. Reliab. Eng.
Syst. Safe., 100:75–83.
9. Chodas, P. W. (1993). Estimating the impact probability of a minor planet with the Earth. Bulletin of the
American Astronomical Society, 25:1236.
10. Chodas, P. W. and Yeomans, D. K. (1999). Could asteroid 1997 XF11 collide with Earth after 2028? In
AAS/Division of Dynamical Astronomy Meeting, volume 31 of Bulletin of the American Astronomical Society, page 1227.
11. Dagum, P., Karp, R., Luby, M., and Ross, S. (2000).
An optimal algorithm for Monte Carlo estimation.
SIAM J. Comput., 29(5):1484–1496.
12. Di Lizia, P. (2008). Robust Space Trajectory
and Space System Design using Differential Algebra.
phdthesis, Politecnico di Milano.
13. Di Lizia, P., Armellin, R., and Lavagna, M. (2008).
Application of high order expansion of two-point
boundary value problems to astrodynamics. Celestial Mechanics and Dynamical Astronomy, 102(4):355–
375.
14. Farnocchia, D., Chesley, S. R., Milani, A., Gronchi,
G. F., and Chodas, P. W. (2015). Orbits, long-term
predicitions, and impact monitoring. In Michel, P., DeMeo, F. E., and Bottke, W. F., editors, Asteroids IV,
pages 815–834. University of Arizona, 1 edition.
15. Klinkrad, H., Alarcon, J. R., and Sanchez, N. (2005).
Collision avoidance for operational ESA satellites. In
Proc. of the 4th European Conference on Space Debris.
16. Kolchin, E. R. (1973). Differential Algebra and Algebraic Groups. Academic Press.
17. Losacco, M., Di Lizia, P., Armellin, R., and Wittig, A. (2018). A differential algebra-based importance
sampling method for impact probability computation
on Earth resonant returns of Near Earth Objects. MNRAS, 479(4):5474–5490.
18. Milani, A. (1999). The asteroid identification problem: I. recovery of lost asteroids. Icarus, 137(2):269–
292.
19. Milani, A., Chesley, S. R., Chodas, P. W., and Valsecchi, G. B. (2002). Asteroid close approaches: Analysis
and potential impact detection. In Bottke Jr., W. F.,
Cellino, A., Paolicchi, P., and Binzel, R. P., editors,
Asteroids III, pages 55–69. The University of Arizona
Press, 1 edition.
20. Morselli, A., Armellin, R., Di Lizia, P., and Bernelli
Zazzera, F. (2014). A high order method for orbital conjunctions analysis: Monte Carlo collision probability
computation. Advances in Space Research, 55(1):311–
33.
21. Patera, R. P. (2001). General method for calculating satellite collision probability. Journal of Guidance,
Control, and Dynamics, 24(4):716–722.
22. Rasotto, M., Morselli, A., Wittig, A., Massari, M.,
Di Lizia, P., Armellin, R., Valles, C. Y., and Ortega, G.
(2016). Differential algebra space toolbox for nonlinear uncertainty propagation in space dynamics. In 6th
International Conference on Astrodynamics Tools and
Techniques (ICATT).
23. Risch, R. H. (1969). The problem of integration in finite terms. Transactions of the American Mathematical
Society, 139:167–189.
24. Risch, R. H. (1970). The solution of the problem
of integration in finite terms. Bulletin of the American
Mathematical Society, 76(3):605–608.
25. Ritt, J. F. (1932). Differential Equations From the
Algebraic Standpoint. American Mathematical Society.
26. Ritt, J. F. (1948). Integration in finite terms: Liouville’s theory of elementary methods. Columbia University Press.
27. Rosenblatt, M. (1952). Remarks on a multivariate
transformation. Ann. Math. Statist., 23(3):470–472.
28. Wittig, A., Di Lizia, P., Armellin, R., Makino, K.,
Bernelli Zazzera, F., and Berz, M. (2015). Propagation
of large uncertainty sets in orbital dynamics by automatic domain splitting. Celestial Mechanics and Dynamical Astronomy, 122(3):239–261.
29. Yeomans, D. K. and Chodas, P. W. (1994). Predicting
close approaches of asteroids and comets to Earth. In
Gehrels, T., Matthews, M. S., and Schumann, A. M.,
editors, Hazards Due to Comets and Asteroids, page
241.
30. Zio, E. (2013). The Monte Carlo Simulation Method
for System Reliability and Risk Analysis. Springer Series in Reliability Engineering. Springer, 1 edition.
31. Zio, E. and Pedroni, N. (2009). Subset Simulation
and Line Sampling for advanced Monte Carlo reliability analysis. In Proceedings of the European Safety and
RELiability (ESREL) 2009 Conference, pages 687–694.
32. Zuev, K. M., Beck, J. L., Au, S.-K., and Katafygiotis, L. S. (2012). Bayesian post-processor and other enhancements of Subsets Simulation for estimating failure probabilities in high dimensions. Computers and
Structures, 92-92:283–296.
</ADVANCED MONTE CARLO SAMPLING TECHNIQUES FOR ORBITAL
CONJUNCTIONS ANALYSIS AND NEAR EARTH OBJECTS IMPACT PROBABILITY
COMPUTATION>