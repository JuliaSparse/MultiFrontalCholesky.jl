using MultiFrontalCholesky, Graphs

function sparseToAdj(M)
    if M.m != M.n
        error("Expected a symmetric matrix")
    end
    n = M.n
    g = simple_adjlist(n,is_directed=false)
    for j=1:n
        for k=M.colptr[j]:M.colptr[j+1]-1
            i = M.rowval[k]
            if i > j
                add_edge!(g,i,j)
            end
        end
    end
    g
end

function laplacian2d(nx,ny)
    M = speye(nx*ny)
    for x=1:nx
        for y=1:ny
            s = x+(y-1)*nx
            M[s,s] = 4
            if x > 1 
                M[s,s-1] = -1
            end
            if x < nx 
                M[s,s+1] = -1 
            end
            if y > 1
                M[s,s-nx] = -1 
            end
            if y < ny 
                M[s,s+nx] = -1
            end
        end
    end
    M
end

(nx,ny) = (100,100)
leafSize = 64;

A = laplacian2d(nx,ny)

g = sparseToAdj(A);
b = randn(nx*ny,1) ;
p = collect(1:size(A,1))

# Bissect
root = Supernode(g,p)
Aperm = A[p,p]
# Build front
front = Front{Float64}(Aperm,root)
# Factorize
MultiFrontalCholesky.Cholesky!(front)

# Check LDLt front
# LD = MultiFrontalCholesky.Unpack(front,root)
# D = diagm(diag(LD))
# L = tril(LD,-1)+speye(size(D,1),size(D,2))
# println("Facto quality : ", norm(Aperm-L*D*L',Inf)/norm(A,Inf))

# Solve Ax = b
sol = MultiFrontalCholesky.Solve(front,root,p,b)
println("Sol quality : ", norm(A*sol-b)/norm(b)) 
#sol2 = A\b 
#println("Done :-)")




