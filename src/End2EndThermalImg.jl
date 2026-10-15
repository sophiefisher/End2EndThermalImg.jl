module End2EndThermalImg

    # surrogate.jl
    export get_transmission
    export compute_and_save_surrogate_transmission_matrix

    # prepare.jl
    export PhysicsHyperParams, ImagingHyperParams, OptimizeHyperParams, ReconstructionHyperParams, JobHyperParams
    export UniformlyRandomObject, UniformlyRandomTFixedDepthObject, GaussianObject
    export FixedNoise, RelativeNoise, apply_noise
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
    using SlurmClusterManager
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
    using BenchmarkTools
    using LinearMaps

    const c = 299792458
    const h = 6.62607015e-34
    const kB = 1.380649e-23
    const grcwa = Ref{Py}()
    const numpy = Ref{Py}()
    
    function __init__()
        grcwa[] = pyimport("grcwa")
        numpy[] = pyimport("numpy")
    end

    include("prepare.jl")
    include("distributed.jl")
    include("forward.jl")
    include("forward_distributed.jl")
    include("backward.jl")
    include("backward_distributed.jl")
    include("optimize.jl")
    include("process.jl")
    include("surrogate.jl")
    include("inversedesign.jl")
    include("diagnostics.jl")

end
