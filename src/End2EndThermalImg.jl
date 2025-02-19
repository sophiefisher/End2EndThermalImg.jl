module End2EndThermalImg

    # prepare.jl
    export PhysicsHyperParams, ImagingHyperParams, OptimizeHyperParams, ReconstructionHyperParams, JobHyperParams

    #surrogate.jl
    export get_transmission

    using CSV
    using DataFrames
    using Interpolations
    using FastChebInterp
    using Parameters
    using PythonCall
    using Dates
    using FFTW
    using Memoization

    const grcwa = Ref{Py}()
    const numpy = Ref{Py}()
    function __init__()
        grcwa[] = pyimport("grcwa")
        numpy[] = pyimport("numpy")
    end

    export grcwa
    export numpy

    include("prepare.jl")
    include("forward.jl")
    include("backward.jl")
    include("optimize.jl")
    include("process.jl")
    include("surrogate.jl")

end
