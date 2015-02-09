type Supernode
    start::Int         # first index of this supernode
    size::Int          # number of indices in this supernode
    struct::Array{Int} # factored structure
    
    children::Array{Supernode}
    has_parent::Bool
    parent::Supernode
    
    function Supernode(g,indices,start=1,parent=nothing)
        this = new()
        this.has_parent = (parent != nothing)
        if parent != nothing
            this.parent = parent
        end
        
        struct_set = IntSet()
        
        # Check if this diagonal block is small enough to treat as dense
        size = length(indices)
        if size <= leafSize
            this.start = start
            this.size = size
            # Push in the original structure
            for s=1:this.size
                for jSub in g.adjlist[s]
                    if jSub > this.size
                        push!(struct_set,jSub+(start-1))
                    end
                end
            end
            this.children = Array(Supernode,0)
        else
            ind_copy = copy(indices)
            
            gL, gR, p, p_inv = bisect(g)
            sizeL = length(gL.vertices)
            sizeR = length(gR.vertices)
            this.start = start + sizeL + sizeR
            this.size = size - (sizeL+sizeR)
            
            indL = sub(indices,1:sizeL)
            indR = sub(indices,sizeL+1:sizeL+sizeR)
            indS = sub(indices,sizeL+sizeR+1:size)
            for k=1:this.size
                indS[k] = ind_copy[p[sizeL+sizeR+k]]
            end

            # Push in the original structure
            for k=1:this.size
                s = p[sizeL+sizeR+k]
                for jSub in g.adjlist[s]
                    if jSub > size
                        push!(struct_set,jSub+(start-1))
                    end
                end
            end
            
            this.children = Array(Supernode,2)
            for k=1:sizeL
                indL[k] = ind_copy[p[k]]
            end
            for k=1:sizeR
                indR[k] = ind_copy[p[k+sizeL]]
            end
            this.children[1] = Supernode(gL,indL,start,      this)
            this.children[2] = Supernode(gR,indR,start+sizeL,this)
        end

        # Perform the symbolic factorization now that the children are formed
        # (recall that the original structure has already been inserted)
        for c in this.children
            for s in c.struct
                if s >= start + size
                    push!(struct_set,s)
                end
            end
        end
        this.struct = collect(struct_set)
        this
    end
end
