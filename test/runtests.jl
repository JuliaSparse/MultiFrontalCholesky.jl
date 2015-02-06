using MultiFrontalCholesky,Compat,Graphs
using Base.Test

function appendel(I,J,V,i,j,v)
    push!(I,i)
    push!(J,j)
    push!(V,v)
end

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

const A = laplacian2d(40,40);

p = [1:size(A,1)]
const root = Supernode(simple_adjlist(A),p);




