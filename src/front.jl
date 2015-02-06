type Front{F}
    L::Array{F,2}
    BR::Array{F,2}
    
    children::Array{Front{F}}
    child_maps::Array{Any,1}
    has_parent::Bool
    parent::Front{F}
    
    function Front(APerm::SparseMatrixCSC{F},snode::Supernode,parent=nothing)
        this = new()
        this.has_parent = (parent != nothing)
        if this.has_parent
            this.parent = parent
        end
        
        # Initialize this node of the matrix
        structSize = length(snode.struct)
        this.L = zeros(F,snode.size+structSize,snode.size)
        for jSub=1:snode.size
            j = jSub + (snode.start-1)
            for k=APerm.colptr[j]:APerm.colptr[j+1]-1
                i = APerm.rowval[k]
                if i >= snode.start && i < snode.start+snode.size
                    iSub = i - (snode.start-1)
                    this.L[iSub,jSub] = APerm.nzval[k]
                    this.L[jSub,iSub] = conj(APerm.nzval[k])
                elseif i >= snode.start+snode.size
                    structRan = searchsorted(snode.struct,i)
                    if length(structRan) == 1
                        iSub = structRan[1] + snode.size
                        this.L[iSub,jSub] = APerm.nzval[k]
                    else
                        error("Found ",length(structRan)," instances of ($i,$j) in struct(",
                        snode.start,":",snode.start+snode.size-1,")=",snode.struct)
                    end
                end
            end
        end
        this.BR = zeros(F,structSize,structSize)
        
        # Handle the children
        num_children = length(snode.children)
        this.children = Array(Front{F},num_children)
        this.child_maps = Array{Any,1}[]
        for c=1:num_children
            this.children[c] = Front{F}(APerm,snode.children[c],this)
            push!(this.child_maps,Array(Int,length(snode.children[c].struct)))
            for k=1:length(snode.children[c].struct)
                i = snode.children[c].struct[k]
                if i < snode.start+snode.size
                    this.child_maps[c][k] = i - (snode.start-1)
                else
                    loc = searchsorted(snode.struct,i)
                    if length(loc) == 1
                        this.child_maps[c][k] = loc[1] + snode.size
                    else
                        error("expected a single index but found ",loc)
                    end
                end
            end
        end

        this
    end 
end

function Unpack!{F}(front::Front{F},snode::Supernode,L::SparseMatrixCSC{F})
    num_children = length(front.children)
    for c=1:num_children
        Unpack!(front.children[c],snode.children[c],L) 
    end
    
    start = snode.start
    totalSize, nodeSize = size(front.L)
    diagInd = start:start+nodeSize-1
    struct = snode.struct
    structSize = length(struct)

    if nodeSize > 0
        L[diagInd,diagInd] = copy(front.L[1:nodeSize,:])
        if structSize > 0
            L[struct,diagInd] = copy(front.L[nodeSize+1:end,:])
        end
    end
end

function Unpack{F}(front::Front{F},snode::Supernode)
    n = snode.start + snode.size-1
    L = speye(n,n)
    Unpack!(front,snode,L)
    L
end

function Cholesky!{F}(front::Front{F})
    # Recurse on the children
    num_children = length(front.children)
    for c=1:num_children
        Cholesky!(front.children[c])
    end
    
    m,nodeSize = size(front.L)
    structSize = m - nodeSize
    FTL = sub(front.L,1:nodeSize,1:nodeSize)
    FBL = sub(front.L,nodeSize+1:m,1:nodeSize)
    FBR = front.BR
    
    # Perform the extend-adds
    for c=1:num_children
        childStructSize = length(front.child_maps[c]);
        for jChild=1:childStructSize
            jSub = front.child_maps[c][jChild];
            for iChild=jChild:childStructSize
                iSub = front.child_maps[c][iChild];
                value = front.children[c].BR[iChild,jChild];
                if iSub <= nodeSize
                    FTL[iSub,jSub] += value;
                elseif jSub <= nodeSize
                    FBL[iSub-nodeSize,jSub] += value;
                else
                    FBR[iSub-nodeSize,jSub-nodeSize] += value;
                end
            end
        end
        # TODO: Clear front.children[c].BR
    end
    
    cholfact!(FTL,:L)
    BLAS.trsm!('R','L','T','N',1.,FTL,FBL)
    BLAS.syrk!('L','N',-1.,FBL,1.,FBR)
end
