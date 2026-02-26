module End2EndThermalImg

    # surrogate.jl
    export get_transmission
    export compute_and_save_surrogate_transmission_matrix

    # prepare.jl
    export PhysicsHyperParams, ImagingHyperParams, OptimizeHyperParams, ReconstructionHyperParams, JobHyperParams
    export UniformlyRandomObject, UniformlyRandomTFixedDepthObject, GaussianObject
    export flatten_object, get_wavcen

    # forward.jl

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
    using PythonCall
    using LaTeXStrings
    using Distributed
    using QuadratureRules
    using InteractiveUtils
    using Zygote
    using ChainRulesCore
    using NLopt
    using JLD2
    using ColorSchemes
    using Statistics
    using Random
    using LinearAlgebra
    using Base.Threads

    const c = 299792458
    const ħ = 6.62607015e-34
    const kB = 1.380649e-23
    const grcwa = Ref{Py}()
    const numpy = Ref{Py}()
    
    function __init__()
        grcwa[] = pyimport("grcwa")
        numpy[] = pyimport("numpy")
    end

    include("prepare.jl")
    include("forward.jl")
    include("backward.jl")
    include("optimize.jl")
    include("process.jl")
    include("surrogate.jl")
    include("inversedesign.jl")
    include("diagnostics.jl")

end
