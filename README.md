# Using-SOSTOOLS-To-Discover-Lyapunov-Functions
Project Summary

This project investigates how Sum-of-Squares (SOS) optimisation methods can be used to automatically discover Lyapunov functions — mathematical tools used to formally prove the stability of dynamical systems.

Lyapunov functions are important because they provide a stability guarantee without needing to solve a system's differential equations directly. Finding them, however, is notoriously difficult. SOS methods offer a systematic approach by converting the problem into a tractable optimisation program, though they come with significant limitations in scalability.

The project was motivated by recent AI breakthroughs in Lyapunov function discovery. SOS methods serve as a rigorous mathematical benchmark against which those emerging AI techniques can be evaluated.

What I Did:

**Benchmarking SOS methods** — I applied an SOS solver (SOSTOOLS in MATLAB) to five dynamical systems of increasing complexity, and to two larger synthetic datasets (200 and 10,000 systems respectively). Candidate Lyapunov functions were independently verified using an SMT solver (Z3 in Python).

**Power system stability analysis** — I applied SOS methods to a real power system model (a two-machine vs infinite-bus system), focusing on estimating the Region of Attraction (ROA): the set of initial conditions from which the system is guaranteed to return to equilibrium. This is increasingly relevant as renewable energy integration reduces natural grid inertia and makes power systems harder to stabilise.

**Comparing ROA algorithms** — I implemented and compared two published algorithms for expanding the ROA estimate, evaluating the trade-off between accuracy and computational cost.

**Scaling to a higher-dimensional system** — I applied the more efficient algorithm to a 10th-order nine-bus AC/DC power system with wind generation, pushing the limits of what SOS methods can handle.

**Key Findings**
SOS methods successfully found and verified Lyapunov functions for 85% of systems in a purpose-built dataset, but only 0.03% of a more general randomly generated dataset — highlighting how sensitive they are to problem structure.
The methods are fundamentally restricted to polynomial systems and cannot discover non-polynomial solutions, even for simple cases.
Computational cost scales poorly with system dimension, making these methods impractical for large real-world systems without approximation.
A more efficient ROA algorithm (Liu et al., 2024) completed in under 3 minutes versus 10–20 minutes for the standard approach, at the cost of a more conservative estimate — a worthwhile trade-off for higher-dimensional systems.
These limitations make the case for AI-driven Lyapunov discovery, while SOS methods remain valuable as a tool for mathematically rigorous verification of AI-generated candidates.

ROA Results: 
<img width="769" height="520" alt="image" src="https://github.com/user-attachments/assets/86dc087b-4b91-4a7d-8f59-dab83d509816" />
<img width="782" height="454" alt="image" src="https://github.com/user-attachments/assets/e40821d5-782c-448f-a6c6-ee7581fd9831" />

 

**Tools & Skills**
MATLAB — SOS optimisation, semidefinite programming, power system modelling
Python — SMT-based formal verification (Z3)
Mathematical methods — Lyapunov stability theory, polynomial optimisation, Region of Attraction estimation
Research skills — literature review, benchmarking, reproducing and comparing published algorithms
Acknowledgements

Carried out as part of an EPSRC-funded summer internship at the University of Liverpool. Builds on work by Alfarano et al. (2024), Anghel et al. (2013), and Liu et al. (2024).
