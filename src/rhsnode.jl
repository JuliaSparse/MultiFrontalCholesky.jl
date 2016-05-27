type RHSNode{F}
    T::Array{F,2}
    B::Array{F,2}
    
    children::Array{RHSNode{F}}
    has_parent::Bool
    parent::RHSNode{F}
    
    function RHSNode(X,snode::Supernode,parent=nothing)
        this = new()
        this.has_parent = (parent != nothing)
        if this.has_parent
            this.parent = parent
        end
        
        n,numRHS = size(X)
        this.T = copy(X[snode.start:snode.start+snode.size-1,:])
        this.B = zeros(F,0,0)
        
        num_children = length(snode.children)
        this.children = Array(RHSNode{F},num_children)
        for c=1:num_children
            this.children[c] = RHSNode{F}(X,snode.children[c],this)
        end
        
        this
    end
end

function Unpack!{F}(XNode::RHSNode{F},snode::Supernode,X::StridedMatrix{F})
    ind = snode.start:snode.start+snode.size-1
    X[ind,:] = copy(XNode.T)
    
    num_children = length(XNode.children)
    for c=1:num_children
        Unpack!(XNode.children[c],snode.children[c],X) 
    end
end

function Unpack{F}(XNode::RHSNode{F},snode::Supernode)
    n = snode.start+snode.size-1
    rootSize, numRHS = size(XNode.T)
    X = zeros(F,n,numRHS)
    Unpack!(XNode,snode,X)
    X
end

function ForwardSolve!{F}(front::Front{F},snode::Supernode,XNode::RHSNode{F})
    # Recurse on the children
    num_children = length(front.children)
    for c=1:num_children
        ForwardSolve!(front.children[c],snode.children[c],XNode.children[c])
    end
    
    # Accumulate the updates from the children
    totalSize, nodeSize = size(front.L)
    nodeSize, numRHS = size(XNode.T)
    structSize = totalSize - nodeSize
    XNode.B = zeros(structSize,numRHS)
    for c=1:num_children
        childStructSize = length(front.child_maps[c]);
        for iChild=1:childStructSize
            iSub = front.child_maps[c][iChild];
            for j=1:numRHS
                value = XNode.children[c].B[iChild,j];
                if iSub <= nodeSize
                    XNode.T[iSub,j] += value;
                else
                    XNode.B[iSub-nodeSize,j] += value;
                end
            end
        end
        # TODO: Clear XNode.children[c].B
    end
    
    LT = sub(front.L,1:nodeSize,1:nodeSize)
    LB = sub(front.L,nodeSize+1:totalSize,1:nodeSize)
    
    # BLAS.trsm!('L','L','N','N',1.,LT,XNode.T)
    # BLAS.gemm!('N','N',-1.,LB,XNode.T,1.,XNode.B)

    BLAS.trsm!('L','L','N','U',1.,LT,XNode.T) # 'U' indicates ones on the diagonal (-> skip diag filled with D)
    BLAS.gemm!('N','N',-1.,LB,XNode.T,1.,XNode.B)

end

function BackwardSolve!{F}(front::Front{F},snode::Supernode,XNode::RHSNode{F})
    totalSize, nodeSize = size(front.L)
    nodeSize, numRHS = size(XNode.T)
    structSize = totalSize - nodeSize
    XNode.B = zeros(structSize,numRHS)
    
    # Pull updates down from the parent
    if front.has_parent
        parentNodeSize = snode.parent.size
        num_siblings = length(front.parent.children)
        
        # Determine which sibling we are
        whichSibling = -1
        for s=1:num_siblings
            if front === front.parent.children[s]
                whichSibling = s
            end
        end
        if whichSibling == -1
            error("This front was not a child of the parent")
        end
        
        for iSub=1:structSize
            iParent = front.parent.child_maps[whichSibling][iSub];
            for j=1:numRHS
                value = XNode.B[iSub,j];
                if iParent <= parentNodeSize
                    XNode.B[iSub,j] += XNode.parent.T[iParent,j]
                else
                    XNode.B[iSub,j] += XNode.parent.B[iParent-parentNodeSize,j]
                end
            end
        end
        # TODO: Clear XNode.parent.B
    end
    
    LT = sub(front.L,1:nodeSize,1:nodeSize)
    LB = sub(front.L,nodeSize+1:totalSize,1:nodeSize)
    
    # BLAS.gemm!('T','N',-1.,LB,XNode.B,1.,XNode.T)
    # BLAS.trsm!('L','L','T','N',1.,LT,XNode.T)
    BLAS.gemm!('T','N',-1.,LB,XNode.B,1.,XNode.T)
    BLAS.trsm!('L','L','T','U',1.,LT,XNode.T) # 'U' indicates ones on the diagonal -> skip
    
    # Recurse on the children
    num_children = length(front.children)
    for c=1:num_children
        BackwardSolve!(front.children[c],snode.children[c],XNode.children[c])
    end
end

function DiagSolve!{F}(front::Front{F},snode::Supernode,XNode::RHSNode{F})

    nodeSize, numRHS = size(XNode.T)

    num_children = length(front.children)
    for c=1:num_children
        DiagSolve!(front.children[c],snode.children[c],XNode.children[c])
    end

    LT = sub(front.L,1:nodeSize,1:nodeSize)
    for row=1:nodeSize
        XNode.T[row,:] = XNode.T[row,:] / LT[row,row] 
    end
end

function Solve{F}(front::Front{F},root::Supernode,p::Array{Int},B::StridedMatrix{F})
    # P A P^T = L L^H and A x = b imply
    #
    #     x = inv(P) (L^H \ (L \ (P b)))
    #
    XNodal = RHSNode{F}(B[p,:],root)
    ForwardSolve!(front,root,XNodal)
    DiagSolve!(front,root,XNodal)
    BackwardSolve!(front,root,XNodal)
    XPerm = Unpack(XNodal,root)
    X = XPerm[invperm(p),:]
    X
end
