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
    using IterativeSolvers
    using Optimisers

    const c = 299792458
    const h = 6.62607015e-34
    const kB = 1.380649e-23
    # Lower bound on T (Kelvin) in reconstruct_Tmap_and_zmap[_distributed]'s NLopt optimization.
    # eps() (~2.22e-16) overflows dB_dT/d2B_dT2's exp(A*freq/T) term (Inf/Inf = NaN) at the
    # wavelengths this package models; confirmed numerically that 3.0 K is safely finite there
    # (though the gradient itself has already underflowed to exactly 0 by that point -- a pixel
    # pinned at this bound still isn't a trustworthy free variable, see T_LOWER_BOUND's use in
    # train_outer_optimization_distributed, src/optimize.jl).
    const T_LOWER_BOUND = 3.0
    const grcwa = Ref{Py}()
    const numpy = Ref{Py}()
    
    function __init__()
        grcwa[] = pyimport("grcwa")
        numpy[] = pyimport("numpy")
        # Scales to however many threads this process itself was started with -- never hardcoded,
        # so this is a no-op (1 thread) on a single-threaded process and uses whatever's actually
        # available otherwise, on the cluster or a laptop alike. FFTW otherwise defaults to 1
        # thread regardless of Threads.nthreads().
        FFTW.set_num_threads(Threads.nthreads())
    end

    include("prepare.jl")
    include("distributed.jl")
    include("forward.jl")
    include("forward_distributed.jl")
    include("backward.jl")
    include("backward_distributed.jl")
    include("optimize.jl")
    include("tune_regularization.jl")
    include("process.jl")
    include("surrogate.jl")
    include("inversedesign.jl")
    include("diagnostics.jl")

end
