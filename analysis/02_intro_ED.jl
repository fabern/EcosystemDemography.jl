using ModelingToolkit
using MethodOfLines
using DomainSets
using OrdinaryDiffEq
using Plots
using Integrals
using StatsFuns

# Example 1: Constant in time and constant with $x$
# TODO: apparently this is a Sinko-Streifer model: 
# https://esajournals.onlinelibrary.wiley.com/doi/epdf/10.2307/1934533
# https://rsmith.math.ncsu.edu/MA573_F17/Population_Models.pdf

@parameters t, x # x is dbh in cm
@variables n(..)

Dt = Differential(t); Dx = Differential(x)

# Using Table 1 (Lehmann 2015) True
x_0 = 1    # cm
r_i = 100  # m-2 ha-1   # TODO: why is this in per hectare?
g_i = 1    # cm year-1
m_i = 0.02 # year-1

domain =[x ∈ Interval(x_0, 150.0), # cm
         t ∈ Interval(0.0, 200.0)] # years
dx = 1.0 # cm
# domain =[x ∈ Interval(x_0, 40.0), # cm
#         t ∈ Interval(0.0, 90.0)] # years
# dx = 0.3 # cm

ic_bc = [n(0.0,  x)       ~ 0.0,  # exp(-(x - 4.0)^2) + exp(-(x - 6.0)^2),
         n(t,  x_0) * g_i ~ r_i]

eq1 = Dt(n(t, x)) ~ - g_i * Dx(n(t, x)) - m_i*n(t, x)
@named sys1 = PDESystem(eq1, ic_bc, domain, [t, x], [n(t, x)])
#                                           indep., depend. variables

# Discretize in x, but leave t undiscretized to use with DifferentialEquations.jl
discretization = MOLFiniteDifference([x => dx], t, approx_order = 2)

prob1 = discretize(sys1, discretization)
sol = solve(prob1, Tsit5(), saveat = 5.00) # save only every 5 years

# Plot transient solution
grid = get_discrete(sys1, discretization)
anim = @animate for (i_t, t_discrete) in enumerate(sol[t])
    plot(grid[x], sol[n(t, x)][i_t, 1:end], 
         label = "Solution", 
         title = "Year: $t_discrete",
         xlabel = "DBH (cm)", ylabel = "tree density n\n(N per cm)", xlim = (x_0, 150.))
    # include also analytical steady state solution
    plot!(grid[x], r_i / g_i * exp.(-m_i/g_i .* grid[x]),
        label = "Exponential distr. (steady-state)", 
        linestyle = :dash, linecolor = :black)
end
gif(anim, "analysis/01_constant_params_A.gif", fps = 8)





# Example 2: Constant in time and varying with $x$
    # # Some parametrisation to fix growth rate:
    # Dmax = 100   # cm,     Basswood from Bonan Figure 19.5
    # hmax = 30    # m,      Basswood from Bonan Figure 19.5
    # G0   = 188.7 # cm/year Basswood from Bonan Figure 19.5
    # b2   = 2 * (hmax - 1.37) / Dmax # from Bonan eq 19.2
    # b3   = (hmax - 1.37) / Dmax^2   # from Bonan eq 19.2
    # h_f(x) = 1.37 + b2*x - b3*x^2 # height (m) as function of bhd (cm) # Taken from Bonan 2019 eq. 19.1
    # plot(x -> h_f(x), 0, 100, xlabel = "DBH (cm)", ylabel = "h (m)") # Should give Bonan Figure 19.4a (Basswood)

# gi_f(x) = G0 * x * (1 - x*h_f(x)/(Dmax * hmax)) / (2.74 + 3*b2*x - 4*b3*x^2)    # cm year-1 # Taken from Bonan 2019 eq. 19.8
# gi_f(x) = 0.001 * G0     * (1 - x*h_f(x)/(Dmax * hmax)) / (2.74 + 3*b2*x - 4*b3*x^2)    # cm year-1 # Taken from Bonan 2019 eq. 19.8
# gi_f(x) = x>20.0 ? 2.00 : 1.00     # cm year-1   # # increase growth rate stepwise above 20 cm 
# @register_symbolic gi_f(x) # special case for step (from https://www.youtube.com/watch?v=8gLhaWRYvfQ)
gi_f(x) = g_i + 0.5 * logistic(x - 30) # cm year-1      # # smoothly increase growth rate above 20 cm
# mi_f(x) = x>20.0 ? m_i : 0.00  # year-1   # # only start dying above 20 cm 
# @register_symbolic mi_f(x) # special case for step (from https://www.youtube.com/watch?v=8gLhaWRYvfQ)
mi_f(x) = m_i + 0 * logistic(x - 20) # year-1   # # only start dying above 20 cm (smoothly)

# Show the varying growth rates and mortality rates
# plot(x -> gi_f(x), 0, 100, xlabel = "DBH (cm)", ylabel = "dDBH/dt (cm/year)") # Should give Bonan Figure 19.4a (Basswood)
# plot(x -> mi_f(x), 0, 100, xlabel = "DBH (cm)", ylabel = "mortality rate (1/year)") # Should give Bonan Figure 19.4a (Basswood)


eq2 = Dt(n(t, x)) ~ - Dx(gi_f(x) * n(t, x)) - mi_f(x) * n(t, x)
@named sys2 = PDESystem(eq2, ic_bc, domain, [t, x], [n(t, x)])
prob2 = discretize(sys2, discretization)
sol2 = solve(prob2, Tsit5(), saveat = 5.00) # save only every 5 years


# prepare analytical solution:
my_int(x_upper) = solve(IntegralProblem((y, p) -> mi_f(y)/gi_f(y), (x_0, x_upper)), 
                            QuadGKJL()).u
my_int.(grid[x])
# plot(1.0:40, my_int.(1.0:40)) # TODO: when 

# Plot
anim2 = @animate for (i_t, t_discrete) in enumerate(sol2[t])
    plot(grid[x], sol2[n(t, x)][i_t, 1:end], 
         label = "Solution", 
         title = "Year: $t_discrete",
         xlabel = "DBH (cm)", ylabel = "tree density n\n(N per cm)", xlim = (x_0, 150.))
    # include also analytical steady state solution
    plot!(grid[x], r_i / g_i * exp.(-m_i/g_i .* grid[x]),
        label = "Exponential distr. (steady-state)", 
        linestyle = :dash, linecolor = :black)
    # additionally the numeric one
    plot!(grid[x], r_i ./ gi_f.(grid[x]) .* exp.(- my_int.(grid[x])),
        label = "Exponential integral (steady-state)", 
        linestyle = :dash, linecolor = :blue)
end
gif(anim2, "analysis/01_constant_params_B.gif", fps = 8)


# Example 2b: Constant in time and varying with $x$
# Test the exponential growth function resulting (Muller-Landau 2006) in Weibull distributions
# TODO...
# gi_f(x) = g_i * x^2.1                    # cm year-1      # # smoothly increase growth rate as power function

    # # # From Muller-Landau 2006 TODO
    # plot!(grid[x], r_i / gi_f(x_0) * grid[x].^(-2.1) .* exp.(-(m_i)/(g_i*(1-2.1)) .* grid[x].^(1-2.1)),
    #     label = "Weibull (Muller-Landau 2006) (steady-state)", 
    #     linestyle = :dash, linecolor = :red)




# Example 3: Constant in time and varying with $x$

# Derive BA diagnostic:
β = 30 # m2 ha-1
xmax = maximum(grid[x])
Ix = Integral(x in DomainSets.ClosedInterval(xmax, x)) # basically cumulative sum from x to xmax
@variables u(..) BA(..) gi_fBA(..)
eq3 = [
    BA(t, x) ~ -Ix(n(t, x)), # Note wrapping the argument to the derivative with an auxiliary variable BA
    gi_fBA(t, x) ~ g_i * max(0, 1 - BA(t, x)/β), # Note wrapping the argument to the derivative with an auxiliary variable BA
    # Dt(n(t, x)) + 2 * n(t, x) + 5 * Dx(BA(t, x)) ~ 1
    Dt(n(t, x)) ~ - Dx(gi_fBA(t, x) * n(t, x)) - mi_f(x) * n(t, x)
]
@named sys3 = PDESystem(eq3, ic_bc, domain, [t, x], [n(t, x)])
prob3 = discretize(sys3, discretization)
sol3 = solve(prob3, Tsit5(), saveat = 5.00) # save only every 5 years


########### INVERSE:
# TODO: since there is an issue with the uppper bound, we can try to use y = -x. So that we would do the lower bound...

@parameters t, y # y is dbh in cm
@variables n(..)
Dt = Differential(t); Dy = Differential(y)
y_0 = -1    # cm
r_i = 100  # m-2 ha-1   # TODO: why is this in per hectare?
g_i = 1    # cm year-1
m_i = 0.02 # year-1
domain =[y ∈ Interval(-150.0, y_0), # cm
         t ∈ Interval(0.0, 200.0)] # years
dy = 1.0 # cm
ic_bc = [n(0.0,  y)       ~ 0.0,
         n(t,  y_0) * g_i ~ r_i]
eq1 = Dt(n(t, y)) ~ g_i * Dy(n(t, y)) - m_i*n(t, y)
@named sys1 = PDESystem(eq1, ic_bc, domain, [t, y], [n(t, y)])
discretization = MOLFiniteDifference([y => dy], t, approx_order = 2)
prob1 = discretize(sys1, discretization)
sol = solve(prob1, Tsit5(), saveat = 5.00) # save only every 5 years

# Plot transient solution
grid = get_discrete(sys1, discretization)
anim = @animate for (i_t, t_discrete) in enumerate(sol[t])
    plot(grid[y], sol[n(t, y)][i_t, 1:end], 
         label = "Solution", 
         title = "Year: $t_discrete",
         xlabel = "DBH (cm)", ylabel = "tree density n\n(N per cm)", xlim = (-150., y_0))
    # include also analytical steady state solution
    plot!(grid[y], r_i / g_i * exp.(m_i/g_i .* grid[y]),
        label = "Exponential distr. (steady-state)", 
        linestyle = :dash, linecolor = :black)
end
gif(anim, "analysis/01_constant_params_A.gif", fps = 8)


using ModelingToolkit
using MethodOfLines
using DomainSets
using OrdinaryDiffEq
using Plots
using Integrals
using StatsFuns
@parameters t, y # y is dbh in cm
@variables n(..)
Dt = Differential(t); Dy = Differential(y)
y_0 = -1    # cm
r_i = 100  # m-2 ha-1   # TODO: why is this in per hectare?
g_i = 1    # cm year-1
m_i = 0.02 # year-1
domain =[y ∈ Interval(-150.0, y_0), # cm
         t ∈ Interval(0.0, 200.0)] # years
dy = 1.0 # cm
ic_bc = [n(0.0,  y)       ~ 0.0,
         n(t,  y_0) * g_i ~ r_i]

discretization = MOLFiniteDifference([y => dy], t, approx_order = 2)
β = 30 # m2 ha-1
ymax = -150 # TODO: same as in grid[]
# gi_f(x) = g_i + 0.5 * logistic(x - 30) # cm year-1      # # smoothly increase growth rate above 20 cm
mi_f(x) = m_i + 0 * logistic(x - 20) # year-1   # # only start dying above 20 cm (smoothly)

Iy = Integral(y in DomainSets.ClosedInterval(ymax, y)) # basically cumulative sum from ymax to y
@variables n(..) BA(..) gi_fBA(..)
# ϵ = 0.001; smooth_larger0(x,ϵ) = log1p(exp(x/ϵ))*ϵ
# smooth_larger0(x,ϵ) = log1p(exp(x/ϵ))*ϵ
b = 0.001; smooth_larger0(arg, b) = exp(asinh(arg/sqrt(b)))*sqrt(b)/2 # https://discourse.julialang.org/t/smooth-approximation-to-max-0-x/109383/14
eq3 = [
    BA(t, y) ~ -Iy(n(t, y)), # Note wrapping the argument to the derivative with an auyiliary variable BA
    # gi_fBA(t, y) ~ g_i * max(0, 1 - BA(t, y)/β), 
    gi_fBA(t, y) ~ g_i * smooth_larger0(1 - BA(t, y)/β, b),
    Dt(n(t, y)) ~ - Dy(gi_fBA(t, y) * n(t, y)) - mi_f(-y) * n(t, y)
]
@named sys3 = PDESystem(eq3, ic_bc, domain, [t, y], [n(t, y), gi_fBA(t, y), BA(t, y)])
prob3 = discretize(sys3, discretization)
sol3 = solve(prob3, Tsit5(), saveat = 5.00) # save only every 5 years

# # Plot transient solution
# grid3 = get_discrete(sys3, discretization)
# anim3 = @animate for (i_t, t_discrete) in enumerate(sol3[t])
#     plot(grid3[y], sol3[n(t, y)][i_t, 1:end], 
#          label = "Solution", 
#          title = "Year: $t_discrete",
#          xlabel = "DBH (cm)", ylabel = "tree density n\n(N per cm)", xlim = (-150., y_0))
#     # include also analytical steady state solution
#     plot!(grid3[y], r_i / g_i * exp.(m_i/g_i .* grid3[y]),
#         label = "Exponential distr. (steady-state)", 
#         linestyle = :dash, linecolor = :black)
# end
# gif(anim3, "analysis/03_BA_A.gif", fps = 8)
###########