# ModelingToolkit and MethodOfLines
# https://www.youtube.com/watch?v=8gLhaWRYvfQ

using ModelingToolkit
using MethodOfLines
using DomainSets

# Example 1: Heat equation with ModelingToolkit and MethodOfLines
@parameters t, x
@variables u(..)

Dt = Differential(t); Dx = Differential(x)
Dxx = Differential(x)^2

α = 1.1
eq = Dt(u(t, x)) ~ α * Dxx(u(t, x))

domain =[x ∈ Interval(0.0, 10.0),
         t ∈ Interval(0.0, 1.0)]

ic_bc = [u(0.0, x) ~ exp(-(x - 4.0)^2) + exp(-(x - 6.0)^2),
         u(t,  0.0) ~ 0.0,
         u(t, 10.0) ~ 0.0]

@named sys = PDESystem(eq, ic_bc, domain, [t, x], [u(t, x)])
#                                         indep., depend. variables

# Discretize in x, but leave t undiscretized to use with DifferentialEquations.jl
dx = 0.1
discretization = MOLFiniteDifference([x => dx], t, approx_order = 2)
prob = discretize(sys, discretization)

using OrdinaryDiffEq
sol = solve(prob, Tsit5(), saveat = 0.05)

# https://sciml.github.io/MethodOfLines.jl/dev/tutorials/brusselator/#Extracting-results
grid = get_discrete(sys, discretization)
using Plots
anim = @animate for (i_t, t_discrete) in enumerate(sol[t])
    plot(grid[x], sol[u(t, x)][i_t, 1:end], 
         ylim = (0., 1.), label = "u", title = "$t_discrete")
end
gif(anim, "analysis/heat_rod.gif", fps = 8)

# Example 2: Heat equation with different initial conditions
step(x) = x>5 ? 1.0 : 0.0
@register_symbolic step(x)
ic_bc2 = [u(0.0, x) ~ step(x),
        u(t,  0.0) ~ 0.0,
        u(t, 10.0) ~ 1.0]
@named sys2 = PDESystem(eq, ic_bc2, domain, [t, x], [u(t, x)])
prob2 = discretize(sys2, discretization)
sol2 = solve(prob2, Tsit5(), saveat = 0.05)

grid2 = get_discrete(sys2, discretization)
anim2 = @animate for (i_t, t_discrete) in enumerate(sol2[t])
    plot(grid2[x], sol2[u(t, x)][i_t, 1:end], 
         ylim = (0., 1.), label = "u", title = "$t_discrete")
end
gif(anim2, "analysis/heat_rod2.gif", fps = 8)


# Example 3: Add advection term
eq3 = Dt(u(t, x)) ~ α * Dxx(u(t, x)) + Dx(u(t, x)) # with advection term

# Example 4: Add time dependent source term
eq4 = Dt(u(t, x)) ~ α * Dxx(u(t, x)) + Dx(u(t, x)) + sin(t)

# Example 5: Test other boundary conditions (e.g. Neumann)
ic_bc5 = [u(0.0, x) ~ step(x),
        Dx(u(t,  0.0)) ~ 0.0,
        u(t, 10.0) ~ 1.0]
# Example 6: Test other boundary conditions (e.g. periodic)
ic_bc6 = [u(0.0, x) ~ step(x),
          u(t,  0.0) ~ u(t, 10.0)]
domain6 = [x ∈ Interval(0.0, 10.0),
           t ∈ Interval(0.0, 5.0)]
eq6 = Dt(u(t, x)) ~ α * Dxx(u(t, x)) + Dx(u(t, x))
@named sys6 = PDESystem(eq6, ic_bc6, domain6, [t, x], [u(t, x)])
prob6 = discretize(sys6, discretization)
sol6 = solve(prob6, Tsit5(), saveat = 0.05)
grid6 = get_discrete(sys6, discretization)
anim6 = @animate for (i_t, t_discrete) in enumerate(sol6[t])
    plot(grid6[x], sol6[u(t, x)][i_t, 1:end], 
         ylim = (0., 1.), label = "u", title = "$t_discrete")
end
gif(anim6, "analysis/heat_rod6.gif", fps = 8)


# Example 7: Multiple spatial dimensions
@parameters t, x, y
Dy = Differential(y)
Dyy = Differential(y)^2
eq7 = Dt(u(t, x, y)) ~ α * (Dxx(u(t, x, y))  + Dyy(u(t, x, y))) + Dx(u(t, x, y)) - Dy(u(t, x, y))
domain7 =[x ∈ Interval(0.0, 10.0),
          y ∈ Interval(0.0, 10.0),
          t ∈ Interval(0.0, 5.0)]
step7(x, y) = x*y>5 ? 1.0 : 0.0
@register_symbolic step7(x, y)
ic_bc7 = [u(0.0, x,   y  ) ~ step7(x, y),
          u(t,   0.0, y  ) ~ u(t, 10.0, y   ),
          u(t,   x,   0.0) ~ u(t, x,    10.0)]
@named sys7 = PDESystem(eq7, ic_bc7, domain7, [t, x, y], [u(t, x, y)])
dx7 = 0.5
dy7 = 0.5
discretization7 = MOLFiniteDifference([x => dx7, y => dy7], t, approx_order = 2)
prob7 = discretize(sys7, discretization7)
sol7 = solve(prob7, Tsit5(), saveat = 0.05)
grid7 = get_discrete(sys7, discretization7)
anim7 = @animate for (i_t, t_discrete) in enumerate(sol7[t])
    heatmap(#grid7[x], 
            sol7[u(t, x, y)][i_t, 1:end, 1:end], 
            # ylim = (0., 1.), 
            label = "u", title = "$t_discrete")
end
gif(anim7, "analysis/heat_square7.gif", fps = 8)


# Example 8: Using parameters (https://sciml.github.io/MethodOfLines.jl/dev/tutorials/params/#Remake-with-different-parameter-values)
@parameters t, x, y, αpar
eq8 = Dt(u(t, x, y)) ~ αpar * (Dxx(u(t, x, y))  + Dyy(u(t, x, y))) + Dx(u(t, x, y)) - Dy(u(t, x, y))
@named sys8 = PDESystem(eq8, ic_bc7, domain7, [t, x, y], [u(t, x, y)], [αpar, ];
                        defaults = Dict(αpar => 1.1))
prob8 = discretize(sys8, discretization7)
newprob8 = remake(prob8)
using SciMLStructures: replace!, Tunable
replace!(Tunable(), newprob8.p, [2.2])
sol8 = solve(newprob8, Tsit5(), saveat = 0.05)
grid8 = get_discrete(sys8, discretization7)
anim8 = @animate for (i_t, t_discrete) in enumerate(sol8[t])
    heatmap(#grid8[x], 
            sol8[u(t, x, y)][i_t, 1:end, 1:end], 
            # ylim = (0., 1.), 
            label = "u", title = "$t_discrete")
end
gif(anim8, "analysis/heat_square8.gif", fps = 8)

# Example 9: Using a Neural Network (min 18:30 in https://www.youtube.com/watch?v=8gLhaWRYvfQ)
# TODO
# eq9 = Dt(u(t, x, y)) ~ αpar * (Dxx(u(t, x, y))  + Dyy(u(t, x, y))) + Dx(u(t, x, y)) - Dy(u(t, x, y)) + NN([u(t, x, y), x, y], ps)
# loss(ps) = sum(abs2, truesol .- solve(remake(prob, p = ps))
