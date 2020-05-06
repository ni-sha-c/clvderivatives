include("cat.jl")
include("clvs.jl")
using Test
u = rand(100,2)
s = rand(2)
@test isapprox(dCat_fd(u,s), dCat(u,s), rtol=1.e-2)
function test_clvs()
    u = rand(1,2)
	s = [rand(),0.2]
    u_trj = Cat(u, s, 20000)[:,:,1]
    du_trj = dCat(u_trj, s)
    les, clv_trj = clvs(du_trj, 2) 
    return les
end
@test isapprox(test_clvs(), [0.9,-0.9],rtol=0.1)


