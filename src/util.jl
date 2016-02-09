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
    g = SimpleAdjacencyList(false,1:n,nedge,alist)
    return g
end

function partToPerm(sizes,part)
    offsets = cumsum(vcat(0,sizes))
    
    p = zeros(part)
    for j in 1:length(p)
        pv = part[j]
        p[offsets[pv+1] += 1] = j
    end
    
    p
end

function bisect(g)
    n = length(g.vertices)

    nedges = 0 
    alist = [sizehint!(Int[],length(g.adjlist[s])) for s in 1:n] ;
    for s=1:n
        for t in g.adjlist[s]
            if t <= n
                push!(alist[s],t) 
                t > s && (nedges += 1)
            end
        end
    end 
    gPruned = SimpleAdjacencyList(false,1:n,nedges,alist) 

    sizes, part = vertexSep(gPruned)
    p = partToPerm(sizes,part)
    p_inv = invperm(p)
    
    sizeL = Int(sizes[1])
    sizeR = Int(sizes[2])
    
    alistL = [sizehint!(Int[],length(g.adjlist[s])) for s in 1:n] ;
    alistR = [sizehint!(Int[],length(g.adjlist[s])) for s in 1:n] ;
    nedgesL = 0 
    nedgesR = 0 
    for s=1:n
        if part[s] == 0
            sMap = Int(p_inv[s])
            for t in g.adjlist[s]
                if t <= n
                    tMap = Int(p_inv[t])                    
                    push!(alistL[sMap],tMap) 
                else
                    push!(alistL[sMap],t)
                end
                nedgesL += 1
            end
        elseif part[s] == 1
            sMap = Int(p_inv[s]-sizeL)
            for t in g.adjlist[s]
                if t <= n
                    tMap = Int(p_inv[t]-sizeL)
                    push!(alistR[sMap],tMap)
                else
                    push!(alistR[sMap],t-sizeL)
                end
                nedgesR += 1
            end
        end
    end
    gL = SimpleAdjacencyList(true,1:Int(sizeL),nedgesL,alistL) 
    gR = SimpleAdjacencyList(true,1:Int(sizeR),nedgesR,alistR) 

    gL, gR, p, p_inv
end
