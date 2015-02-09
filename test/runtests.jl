using MultiFrontalCholesky,Compat,Graphs
using Base.Test

const A = laplacian2d(40,40);

p = [1:size(A,1)]
const root = Supernode(simple_adjlist(A),p);

@test p[1:5] == [1289,1328,1329,1330,1369]




