module DyadBotComponents

include("blocks.jl") # Blocks parameter definitions

include("DiscreteKalmanFilter.jl")
include("planar_flat.jl")
include("cascade_planar_flat.jl")



# include("planar_multibody.jl")
# include("segway_3d.jl")

end # module DyadBotComponents