A very rough extension to the `DifferentialEquations.jl` ecosystem to solve generic [Landau-Lifshitz-Gilbert equations](https://en.wikipedia.org/wiki/Landau%E2%80%93Lifshitz%E2%80%93Gilbert_equation). 

These problems are defined as
$$
\partial_t \vec{m} = \vec{m} \times \vec{A}(m,t),
$$
where $\vec{m}$ is the magnetization and $\vec{A}$ is a function of the magnetization. 
E.g. for a simple magnetic field along the $z$ direction, $\vec{A}=-\gamma H_0 \hat{z}$ and the spins will start precessing. In general $\vec{A}(m,t)$ can depend on both magnetization and time, making the problem non-linear. 

You have to provide everything yourself using the $\vec{A}(m,t)$ function, including the Gilbert damping like term. Stochastic fluctuations are also supported.