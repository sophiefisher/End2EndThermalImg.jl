module End2EndThermalImg

    # prepare.jl
    export PhysicsHyperParams, ImagingHyperParams, OptimizeHyperParams, ReconstructionHyperParams, JobHyperParams
    export UniformlyRandomObject
    export get_wavcen
    
    #surrogate.jl
    export get_transmission
    export compute_and_save_surrogate_transmission_matrix

    using CSV
    using DataFrames
    using Interpolations
    using FastChebInterp
    using Parameters
    using PythonCall
    using Dates
    using FFTW
    using Memoization
    using FastChebInterp
    using PythonPlot
    using LaTeXStrings

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
