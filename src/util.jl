@doc "Generate a simple adjacency list graph representation from a symmetric or Hermitian sparse matrix"->
function Graphs.simple_adjlist{Tv,Ti}(S::SparseMatrixCSC{Tv,Ti})
    issym(S) || ishermitian(S) || error("S must be symmetric or Hermitian")
    n = size(S,1)
    colpt = S.colptr
    nedge = 0

    alist = [@compat sizehint!(Ti[],colpt[j+1]-colpt[j]) for j in 1:n]
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

@doc """Create a permutation from the partition in `part`.

`sizes` is redundant but already calculated
"""->
function partToPerm(sizes,part)
    offsets = cumsum(vcat(0,sizes))
    
    p = zeros(part)
    for j in 1:length(p)
        pv = part[j]
        p[offsets[pv+1] += 1] = j
    end
    
    p
end

@doc "Bisect a graph using the partition created by vertexSep"->
function bisect(g)
    n = length(g.vertices)
    gPruned = simple_adjlist(n,is_directed=false)
    for s=1:n
        for t in g.adjlist[s]
            if t > s && t <= n
                add_edge!(gPruned,s,t)
            end
        end
    end 
    sizes, part = vertexSep(gPruned)
    p = partToPerm(sizes,part)
    p_inv = invperm(p)
    
    sizeL = Int(sizes[1])
    sizeR = Int(sizes[2])
    
    gL = simple_adjlist(Int(sizeL),is_directed=true)
    gR = simple_adjlist(Int(sizeR),is_directed=true)
    for s=1:n
        if part[s] == 0
            sMap = Int(p_inv[s])
            for t in g.adjlist[s]
                if t <= n
                    tMap = Int(p_inv[t])
                    add_edge!(gL,sMap,tMap)
                else
                    add_edge!(gL,sMap,t)
                end
            end
        elseif part[s] == 1
            sMap = Int(p_inv[s]-sizeL)
            for t in g.adjlist[s]
                if t <= n
                    tMap = Int(p_inv[t]-sizeL)
                    add_edge!(gR,sMap,tMap)
                else
                    add_edge!(gR,sMap,t-sizeL)
                end
            end
        end
    end
    
    gL, gR, p, p_inv
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
