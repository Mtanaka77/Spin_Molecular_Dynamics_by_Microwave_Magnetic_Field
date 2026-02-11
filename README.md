# Spin_molecular_dynamics_by_magnetic_microwaves

## Sintering Experiments and Theory

Microwave sintering by initial experiments (R. Roy et. al., Nature, 1999), Ref.1.  

Theory identifying unpaired 3d electron spins of Fe(3+) and Fe(2+) to increase 
much above the Curie temperature (M. Tanaka et al., Phys. Rev. B, 2009), Ref.2.  


## Simulation Procedures and Our Theory

Electron spins of Fe(3+), Fe(2+) and O(2-) in cubic cells

Magnetic microwaves of giga-Hertz frequency, like 2.5 GHz

Metropolis criterion is used to accept/reject in the MC cycle

Dissipation spin molecular dynamics simulation is executed in the MD cycle.

  > Execution of a few 1,000,000 steps (1 GHz ranges) for $ Delta t= 0.001 $ ps

Microwave magnetic heating by numerical simulation reaches 1,300 Celsius for metal oxide magnetite.
The peak of the time derivative dU_sys/dt corresponds to the Curie temperatre of 858 K (585 C). 
The total energy U_sys and the time derivative dU_sys/dt are shown as the proof of our theory, Ref. 2. 

## Numerical Code, Parameters and Files

1) @spin_SMD11a.f03: Dissipative spin molecular dynamics. An odd number of processors
in the z direction is required for your choice. Read the simulation code and parameters
for details !!

3) param_spinRL11.h: Parameters 

 > number of processors (nodes), total cells of irons and oxygens (odd numbers),
p3m resolution
  
3) SAI111_config.START1: Configuration file

 > physical run time, lattice sizes, number of cells, exchange integrals for Fe and O,
  period of microwave magnetic field, temperature, Curie temprature, etc.
  About 1,000,000 steps are required !

4) Magnetite in a cubic lattice: magnetite8.xyz for initialization

## References

1. R. Roy, D. Agrawal, J. Cheng, and S. Gedevanishvili, Full sintering of powdered-metal bodies
in a magnetic field, Nature, 399, 668 (1999).

2. M. Tanaka, H. Kono, and K. Maruyama, Selective heating mechanism of magnetic
metal oxides by a microwave magnetic field, Phys. Rev. B, 79, 104420 (2009).

