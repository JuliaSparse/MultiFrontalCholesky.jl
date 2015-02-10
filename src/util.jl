@doc "Generate a simple adjacency list graph representation from a symmetric or Hermitian sparse matrix"->
function Graphs.simple_adjlist(S::SparseMatrixCSC)
    issym(S) || ishermitian(S) || error("S must be symmetric or Hermitian")
    n = size(S,1)
    colpt = S.colptr
    nedge = 0

    alist = [@compat sizehint!(Int[],colpt[j+1]-colpt[j]) for j in 1:n]
    for j in 1:n
        for k in colpt[j]:(colpt[j+1]-1)
            i = S.rowval[k]
            if i != j                   # omit self-edges
                push!(alist[j],i)
                i > j && (nedge += 1)
            end
        end
    end
    SimpleAdjacencyList(false,1:n,nedge,alist)
end

immutable PartSet
    left::IntSet
    right::IntSet
    center::IntSet
end

function PartSet(part::Vector)
    l = IntSet(); r = IntSet(); c = IntSet()
    for i in 1:length(part)
        0 <= (p = part[i]) <= 2 || error("part[$i] = $(part[i]) should be in [0,1,2]")
        push!(p == 0 ? l : (p == 1 ? r : c), i)
    end
    PartSet(l,r,c)
end

@doc "bisect the graph g using a vertex separator from Metis"->
function bisect(g::SimpleAdjacencyList)
    ps = PartSet(vertexSep(g)[2])
                                        # create a permutation from ps
    p = vcat(collect(ps.left),collect(ps.right),collect(ps.center))
    invp = invperm(p)                   # also checks that p is a permutation

    nL = length(ps.left)
    nR = length(ps.right)
    nLR = nL + nR
    alistL = [Int[] for _ in 1:nL]
    alistR = [Int[] for _ in 1:nR]
    nedgeL = nedgeR = 0

    for i in 1:nL                       # assemble gL
        ii = p[i]
        for j in g.adjlist[ii]
            nL < (jj = invp[j]) <= nLR && error("ps does not partition g")
            if jj < nL # omit edges between elements of ps.left and ps.center
                push!(alistL[i],jj)
                jj > i && (nedgeL += 1)
            end
        end
    end
    
    for i in 1:nR                       # assemble gR
        ii = p[nL + i]
        for j in g.adjlist[ii]
            (jj = invp[j]) <= nL && error("ps does not partition g")
            if jj < nLR
                jj -= nL
                push!(alistR[i],jj)
                jj > i && (nedgeR += 1)
            end
        end
    end
    SimpleAdjacencyList(false,1:nL,nedgeL,alistL),
    SimpleAdjacencyList(false,1:nR,nedgeR,alistR),p,invp
end

@doc "Push i, j and v onto I, J and V, respectively"->
function appendel(I,J,V,i,j,v)
    push!(I,i)
    push!(J,j)
    push!(V,v)
end

@doc "Generate a 2-dimensional Laplacian matrix for a mesh?"->
function laplacian2d(nx,ny)
    n = nx*ny
    nzest = 5n
    I = @compat sizehint!(Int32[],nzest)
    J = @compat sizehint!(Int32[],nzest)
    V = @compat sizehint!(Float64[],nzest)
    for x in 1:nx
        for y in 1:ny
            s = x + (y-1)*nx
            appendel(I,J,V,s,s,2)
            x > 1 && appendel(I,J,V,s,s-1,-1)
            y > 1 && appendel(I,J,V,s,s-nx,-1)
        end
    end
    A = sparse(I,J,V,n,n)
    A + A'
end
