# Coupling the Phase-Field approach with a Discrete Element Model to investigate Pressure-Solution

TODO add download

## Problem description

Even if most of the materials appear homogeneous, they are heterogeneous at smaller scales.
For geomaterials, this smaller scale is at the level of the microstructure and it consists of grains (that can differ in minerals) and pores.
This microstructure is fundamental as the organization and the elements that compose it determine the multiphysical properties at the scale of the material.
Furthermore, this microstructure is affected by external solicitations and evolves with time.
Then, it becomes essential to predict the microstructure evolution.
A relevant example of this change in porous geomaterials is the pressure-solution phenomenon.
As depicted in [Figure 1], three fundamental chemo-mechanical processes at the microscale are involved: (1) dissolution due to stress concentration at grain contacts, (2) diffusive transport of dissolved mass from the contact to the pore space, and (3) precipitation of the solute on the less stressed surface of the grains.

<a id="fig-configuration"></a>

![Definition of the configuration and description of the microstructure with a phase variable](fig_configuration.png)

***Figure 1:** Definition of the configuration and description of the microstructure with a phase variable.*


### Phase-Field description

Among the methods available in the literature, the phase-field approach emerges as an efficient and accurate tool for this prediction of the microstructure evolution due to multiphysical solicitations. 
This new example is a continuation of the previously modeled [dissolution, diffusion, and precipitation pattern](../case3_dissolution_diffusion_precipitation/case3_dissolution_diffusion_precipitation.md)

The dissolution/precipitation of the solid phase can be described by an Allen-Cahn equation, available in Eq. 1. 
This equation is applied to a phase variable $\eta$ ($=1$ if the point corresponds to the solid phase, $=0$ if the point corresponds to the pore space).

\[
    \frac{\partial\eta}{\partial t}=-L\,\frac{\partial\left(f_{loc}+E_d\right)}{\partial \eta} + L\cdot\kappa\,\nabla^2\eta
    \label{Equation Allen Cahn}
    \tag{1}
\]

In particular, this equation involves a free energy function $f_{loc}$ that describes the material and depends on the mineral. 
This free energy function is destabilized by a tilting energy function $E$, inducing the evolution of the microstructure, as depicted in [Figure 2] and formulated in Eq. 2. 
The additional parameters $L$ and $\kappa$ from Eq. 1 affect the interface width and the kinetics of the phenomenon.

<a id="fig_floc_ed"></a>


![Destabilization of the free energy $f_{loc}$ by the tilting energy $E$](fig_floc_E.png){ width="60%" }

***Figure 2:** Destabilization of the free energy $f_{loc}$ by the tilting energy $E$.*

\[
    f_{loc}+E(c) = W\times\left(\eta^2(1-\eta)^2\right) + ed(c)\times\left( 3\eta^2-2\eta^3\right)
    \label{Equation free energy}
    \tag{2}
\]
where $W$ is the barrier energy, preventing the phase transition, and $ed$ the tilting amplitude, inducing the phase transition. 
The amplitude of the tilting is dictated by the value of a new variable $c$ that describes the concentration of a reactive specimen in the pore fluid.
In the context of this model, the propagation of this solute is described by Eq. 3.

\[
    \frac{\partial c}{\partial t} = -\alpha_\textit{source} \frac{\partial \textit{source}}{\partial t} -\alpha_\textit{product} \frac{\partial \textit{product}}{\partial t} + \kappa_c \nabla^2 c
    \label{Equation Diffusion c}
    \tag{3}
\]

the terms $\frac{\partial c}{\partial t} = \kappa_c \nabla^2 c$ represent the diffusive propagation of the solute, and the terms $\alpha_\eta \frac{\partial \eta}{\partial t}$ ensure the conservation of the mass during the dissolution/precipitation of the solid.

Subsequently, by assuming the chemical reaction $\left(\eta_s \rightleftharpoons c_l\right)$ at a fixed pressure $P$, the chemical quotient of the reaction is $Q=\frac{\{c\}}{\{\eta\}}=c$.
In the same note, the equilibrium constant is $K=c_{eq}(P)$, which is pressure-dependent.
The dissolution of a solid phase occurs for $Q<K$ ($c<c_{eq}$) and the precipitation of a solid phase occurs for $Q>K$ ($c>c_{eq}$).
Furthermore, the kinetics of the reactions depend on the distance to the equilibrium, see Eq. 4.

\[
    E = \chi_{diss/prec}\times a_s\times \left(1 - \frac{c}{c_{eq}\,a_s}\right)\times \left(3\eta^2-2\eta^3\right)
    \label{Equation E}
    \tag{4}
\]
where $\chi_{diss/prec}$ is a constant relative to the global kinetics of the reaction (dissolution or precipitation) and $a_s$ is the solid activity that is defined in Eq. 5.

\[
    a_s = \text{exp}\left(\frac{P\times V_m}{RT}\right)
    \label{Equation as}
    \tag{5}
\]
where P is the pressure at the contact, $V_m$ is the molar volume, $R$ is the gas constant, and $T$ is the temperature.
The solid activity $a_s$ appears to increase in proportion to the stress transmitted at the contact. This modification affects the value of the solute concentration at equilibrium $c_{eq}\, a_s$.
The established chemical equilibrium in the fluid film ($c$ versus $c_{eq}\,a_s$) becomes unverified.
The material dissolves in the contact zone to reach the chemical equilibrium $c=c_{eq}\,a_s$.
A gradient in the solute concentration values appears between the contact zone $c=c_{eq}\,a_s$ and the pore space $c = c_{eq}$.
This gradient gives rise to the diffusion of the solute.
Subsequently, the established chemical equilibrium in the pore space ($c$ versus $c_{eq}$) becomes unverified. The material precipitates in the pore space to reach the chemical equilibrium $c=c_{eq}$.


As illustrated by Figure 1, the dissolution of the solid can remove a pre-existing contact.
Subsequently, the grains, which are loaded, should reorganize to reach a new mechanical steady-state (in the case depicted herein, the particle moves down to retrieve the contact).
It is essential to emphasize that the Phase-Field formulation itself is not able to compute the granular mechanical steady-state and interpolate the stress at the contacts.
Thus, the idea is to couple this approach with another method that has these abilities: the Discrete Element Model.

### Discrete Element Model

The Discrete Element Model has been developed to account for the individual grains and their interactions.
In particular, the momentum balances (see Eq. 6) are solved for each grain, considering the external forces applied to the system $f/M^{external}$ and the transmitted forces interpolated at the contacts $f/M^{contacts}$. 
These contact forces are determined through a penalization method based on an overlap $\delta$ between the grains.
It is essential to specify that this overlap $\delta$ can be based on a characteristic distance or volume of the contact.

$$
\begin{align*}
    m\,\frac{\partial v_i}{\partial t}&=f_i^{external} + f_i^{contacts}(\delta)\nonumber\\
    I\,\frac{\partial w_i}{\partial t}&=M_i^{external} + M_i^{contacts}(\delta)
    \label{Equation Momentum}
    \tag{6}
\end{align*}
$$

where $m$ is the mass of the grain, $v_i$ the particle velocity (with Einstein's notation), $I$ is the moment of inertia of the grain, and $w_i$ is the angular velocity .

In the case depicted in Figure 1, the granular reorganization remains trivial: the grain moves down until it verifies the overlap $\delta$ that corresponds to the external force applied $F$.
In a more complex configuration (which is common in the literature), this reorganization should be estimated with dedicated solvers such as [Yade](https://yade-dem.org/doc/).

### Description of the Phase-Field Discrete Element Method

---

## Model set up

---

## Results



# References

1. N. Moelans, B. Blanpain, P. Wollants (2008), "An introduction to phase-field modeling of microstructure evolution", Computer Coupling of Phase Diagrams and Thermochemistry 32:268–294, DOI: 10.1016/j.calphad.2007.11.003
2. C. O'Sullivan (2011), "Particulate Discrete Element Modeling", CRC Press, DOI: 10.1201/9781482266498
3. A. Guével, H. Rattez, E. Veveakis (2020), "Viscous phase-field modeling for chemo-mechanical microstructural evolution: application to geomaterials and pressure solution", International Journal of Solids and Structures 207:230-249, DOI: 10.1016/j.ijsolstr.2020.09.026
4. A. Sac-Morane, M. Veveakis, H. Rattez (2024), "A Phase-Field Discrete Element Method to study chemo-mechanical coupling in granular materials", Computer Methods in Applied Mechanics and Engineering 424:116900, DOI: 10.1016/j.cma.2024.116900
5. A. Sac-Morane, H. Rattez, M. Veveakis (2025), "Importance of precipitation in the slowdown of creep behavior induced by pressure solution", Journal of Engineering Mechanics 151:04025025, DOI: 10.1061/JENMDT.EMENG-8360