include("cat.jl")
using Test
u = rand(100,2)
s = rand(2)
@test isapprox(dCat_fd(u,s), dCat(u,s), rtol=1.e-6)


