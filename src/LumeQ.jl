module LumeQ

# Include threerd-party dependencies
using LinearAlgebra


# Include all submodules
include("dynamics/dynamics.jl")
include("embedding/embedding.jl")
include("opt/opt.jl")
include("plot/plot.jl")
include("property/property.jl")
include("spins/spins.jl")
include("states/states.jl")
include("transport/transport.jl")
include("utils/utils.jl")


# Make submodules accessible from the main module
using .Dynamics
using .Embedding
using .Opt
using .Plot
using .Property
using .Spins
using .States
using .Transport
using .Utils


# re-export commonly used functions and types


end # module
