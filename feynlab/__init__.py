"""feynlab — a Feynman-diagram laboratory for many-body physics and Hall transport.

Subpackages
-----------
core        Matsubara grids, Fourier transforms, lattices/DOS, Pade continuation
diagrams    enumeration, classification and Monte Carlo evaluation of Feynman diagrams
selfenergy  local (HF, 2OPT, Hubbard-I, atomic) and lattice (2nd order, RPA/FLEX) self-energies
qmc         CT-INT QMC, polaron DiagMC, dense & sector ED, ED-DMFT bath fitting
dmft        Bethe-lattice DMFT with IPT / Hubbard-I / CT-INT solvers
egas        RPA electron gas (Lindhard, plasmons, correlation energy) and BCS
transport   Berry-curvature, Kubo-Streda (ladder + skew) and Boltzmann (side jump, skew)
            Hall responses for charge, spin, orbital and orbital-quadrupole currents
book        companion calculators for Mattuck's 'A Guide to Feynman Diagrams'
"""
__version__ = "0.2.0"
