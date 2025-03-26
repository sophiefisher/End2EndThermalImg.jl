module End2EndThermalImg

    # prepare.jl
    export PhysicsHyperParams, ImagingHyperParams, OptimizeHyperParams, ReconstructionHyperParams, JobHyperParams
    export UniformlyRandomObject
    export get_wavcen
    
    # surrogate.jl
    export get_transmission
    export compute_and_save_surrogate_transmission_matrix

    # python modules
    export grcwa
    export numpy

    using CSV
    using DataFrames
    using Interpolations
    using FastChebInterp
    using Parameters
    using Dates
    using FFTW
    using Memoization
    using FastChebInterp
    using PythonPlot
    using LaTeXStrings

    const c = 299792458
    const ħ = 6.62607015e-34
    const kB = 1.380649e-23

    include("prepare.jl")
    include("forward.jl")
    include("backward.jl")
    include("optimize.jl")
    include("process.jl")
    include("surrogate.jl")

end
