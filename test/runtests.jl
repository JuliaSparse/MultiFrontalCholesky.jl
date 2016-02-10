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

function appendel(I,J,V,i,j,v)
    push!(I,i)
    push!(J,j)
    push!(V,v)
end

function laplacian2d(nx,ny)
    n = nx*ny
    nzest = 5*n
    I = sizehint!(Int32[],nzest)
    J = sizehint!(Int32[],nzest)
    V = sizehint!(Float64[],nzest)
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




