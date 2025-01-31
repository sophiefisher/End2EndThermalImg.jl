module End2EndThermalImg

    # prepare.jl
    export PhysicsHyperParams

    #surrogate.jl
    export get_transmission

    using CSV
    using DataFrames
    using Interpolations
    using FastChebInterp
    using Parameters
    using PythonCall
    using Dates

    include("prepare.jl")
    include("forward.jl")
    include("backward.jl")
    include("optimize.jl")
    include("process.jl")
    include("surrogate.jl")

end
