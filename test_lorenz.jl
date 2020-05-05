include("lorenz.jl")
include("clvs.jl")
using Test
u = rand(100,3)
@test isapprox(dlorenz63(u), dlorenz63_fd(u), atol=1.e-5)
function test_clvs()
	u = rand(1,3)
	s = [10.,28.,8/3.]
	u_trj = lorenz63(u, s, 20000)[:,:,1]
	du_trj = dlorenz63(u_trj, s)
	les, clv_trj = clvs(du_trj, 3) 
	return les
end
@test isapprox(test_clvs()/0.005, [0.9,0.,-14.5],rtol=0.1)
