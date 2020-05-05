using LinearAlgebra
function lorenz63(u::Array{Float64,2},
				  s::Array{Float64,1}=[10.,28.,8/3],
				  n::Int64=1)
	m, d = size(u)
	u_trj = zeros(m,d,n+1) 
	u_trj[:,:,1] = u
	dt = 0.005
	for i = 2:n+1
		x = view(u_trj,:,1,i-1)
		y = view(u_trj,:,2,i-1)
		z = view(u_trj,:,3,i-1)
		x1 = view(u_trj,:,1,i)
		y1 = view(u_trj,:,2,i)
		z1 = view(u_trj,:,3,i)
		@. x1 = x + dt*s[1]*(y - x)
		@. y1 = y + dt*(x*(s[2] - z) - y)
		@. z1 = z + dt*(x*y - s[3]*z)
	end
	return permutedims(u_trj,[3,2,1])
end
function dlorenz63(u::Array{Float64,2}, s::Array{Float64,1}=
			   [10.,28.,8/3],dt::Float64=0.005)
	m, d = size(u)
	dTu = ones(d,d,m) 
	for i=1:m
		dTu[1,1,i] += -dt.*s[1] 
		dTu[1,2,i] = dt.*s[1]
		dTu[1,3,i] = 0.0
		
		dTu[2,1,i] = @. dt*(s[2] - u[i,3]) 
		dTu[2,2,i] += -dt
		dTu[2,3,i] = -dt*u[i,1]
		
		dTu[3,1,i] = dt.*u[i,2]
		dTu[3,2,i] = dt.*u[i,1]
		dTu[3,3,i] += -dt*s[3]
	end
	return dTu
end
function dlorenz63_fd(u::Array{Float64,2}, s::Array{Float64,1}=[10.,28.,8/3],dt::Float64=0.005)
	
	m,d = size(u)
	dTu_fd = ones(d,d,m)
	
	eps = 1.e-4
	u_p = copy(u)
	u_p[:,1] .+= eps
	u_m = copy(u)
	u_m[:,1] .-= eps
	dTx_p = (lorenz63(u_p, s, 1)[2,:,:] - 
			 lorenz63(u_m, s, 1)[2,:,:])/(2.0*eps)
	
	u_p[:,1] .= u[:,1]
	u_p[:,2] .+= eps
	u_m[:,1] .= u[:,1]
	u_m[:,2] .-= eps
	dTy_p = (lorenz63(u_p, s, 1)[2,:,:] - 
			 lorenz63(u_m, s, 1)[2,:,:])/(2.0*eps)

	u_p[:,2] .= u[:,2]
	u_p[:,3] .+= eps
	u_m[:,2] .= u[:,2]
	u_m[:,3] .-= eps
	dTz_p = (lorenz63(u_p, s, 1)[2,:,:] - 
			 lorenz63(u_m, s, 1)[2,:,:])/(2.0*eps)

	dTu_fd = reshape([dTx_p; dTy_p; dTz_p],d,d,m)
	return dTu_fd
end

#=
u = rand(1,3)
s = [10.,28.,8/3.]
@time u_trj = lorenz63(u, s, 20000)[:,:,1]
@time du_trj = dlorenz63(u_trj, s)
@time les, clv_trj = clvs(du_trj, 3) 
=#
