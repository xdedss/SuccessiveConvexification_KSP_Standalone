# SuccessiveConvexification

Implementation of [Successive Convexification for 6-DoF Mars Rocket Powered Landing with Free-Final-Time](https://arxiv.org/abs/1802.03827)

Reqires `python3`, `matplotlib`, `numpy`, `scipy`, `sympy` and `cvxpy 1.0.0`.

Since [cvxpy 1.0.0](https://github.com/cvxgrp/cvxpy/tree/1.0) only supports Mac and Linux, Windows users are reccomended to use the [Windows Subsystem for Linux (WSL) and Ubuntu](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). This is confirmed to work with reasonable speed (this is my development environment).


This branch:
- Erases MATLAB (.m) code for symbolic dynamics code generation
  - Replaced with `sympy` symbolic code generation


- Erases visualization code - for no reason other than I haven't yet used it really :)
  - Replaced with simple matplotlib plots


- Adds a tremendous amount of documentation and pythonic operations
- Mostly conforms to PEP8

- Due to intense hatred of variables ending with `_`, replaced all cvx variables ending with `_` with proper names
   - v instead of _ denotes a cvx Variable
   - parm instead of _ denotes a cvx Parameter


- Replaced a great number of () with [] to clarify definitions of arrays versus function calls

- Replaced all dynamics functions in parameters with the sympy generated ones in dynamics_functions.py

----

# How to Run

Solve dynamics matrices and generate all functions : `python dynamics_generation.py`

Solve successive convex optimization problems : `python main.py`

Plot solution variables : `python ./trajectory/plot.py`


[Short Video Visualization by Sven Niederberger - from the master branch](https://gfycat.com/HideousUniqueEthiopianwolf)

----

Roadmap:

- Lift and Drag Dynamics With Aerodynamic Controls, in the spirit of [this paper by Liu](https://arc.aiaa.org/doi/abs/10.2514/6.2017-1732)

- Automatic scaling of inputs from real-world dimensional scenarios to dimensionless solution inputs, and the reverse for the output

- C++ implementation likely to be forked from EmbersArc / Sven Niederberger's [WIP C++ implementation](https://github.com/EmbersArc/SuccessiveConvexificationCpp)

- Real-time performance evaluation in a realistic spaceflight simulation



History of this project:

- After the paper came out, I immediately started working on a python implementation. I worked for about 2 weeks until I came across [Sven Niederberger's](https://github.com/EmbersArc) implementation on github! I ditched my implementation although it was nearly convergent.
