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