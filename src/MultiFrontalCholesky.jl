module MultiFrontalCholesky
    using Compat,Graphs,Metis

    if VERSION < v"0.4-"
        using Docile                    # for @doc macro
    end

    export Cholesky!,
           Front,
           RHSNode,
           Supernode,
           Solve,
           laplacian2d                  # generate test matrices

    const leafSize = 64
           
    include("util.jl")
    include("supernode.jl")
    include("front.jl")
    include("rhsnode.jl")

end # module
