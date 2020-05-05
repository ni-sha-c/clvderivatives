using LinearAlgebra
"""
	n: number of timesteps
    m: number of initial conditions
    s[0] = abs(lambda), s[1] = alpha
    output size: mxdx(n+1)
"""
function Cat(u::Array{Float64,2}, s::Array{Float64,1}
			  =[0.,0.], n::Int64=1)
	theta(phi) =  2*pi*phi - s[2]
	Psi(phi) =  @. (1/pi)*atan(s[1]*sin(theta(phi))/
						(1. - s[1]*cos(theta(phi))))
	m, d = size(u)
    u_trj = zeros(m,d,n+1)
	u_trj[:,:,1] = u
	
    for i = 1:n
		x = view(u_trj,:,1,i)
		y = view(u_trj,:,2,i)

        psix = Psi(x)
		u_trj[:,1,i+1] = @. (2*x + y + psix) % 1
		u_trj[:,2,i+1] = @. (x + y + psix) % 1
	end
	return permutedims(u_trj,[3,2,1])
end
"""
    n: number of timesteps
    m: number of initial conditions
    s[0] = abs(lambda), s[1] = alpha
    output size: (n+1)xdxm
"""
function inverse_Cat(u::Array{Float64,2}, s::Array{Float64,1}=[0.,0.], n::Int64=1)
	
	Psi2(phi) = @. (1/pi)*atan(s[1]*sin(2*pi*phi-s[2])/
		(1. - s[1]*cos(2*pi*phi - s[2])))
	m, d = size(u)
	u_trj = zeros(m,d,n+1)
	u_trj[:,:,1] = u
	psi2 = zeros(m)
    for i = 1:n
		x = view(u_trj,:,1,i)
		y = view(u_trj,:,2,i)
        psi2 .= Psi2(x .- y)
        u_trj[:,1,i+1] = @. (x - y) % 1
        u_trj[:,2,i+1] = @. (-x + 2*y - psi2) % 1
	end
	return permutedims(u_trj,[3,2,1])
end
function dCat(u::Array{Float64,2}, s::Array{Float64,1}=[0.,0.])
	m, d = size(u)
	theta(phi) =  2*pi*phi - s[2]
    dtheta = 2*pi
	num_t(phi) =  s[1]*sin(theta(phi))
	den_t(phi) = 1. - s[1]*cos(theta(phi))
	t(phi) = num_t.(phi)./den_t.(phi)
	Psi(t) = (1/pi)*atan.(t)
	dnum_t(phi) = s[1]*cos.(theta.(phi))*dtheta
	dden_t(phi) = s[1]*sin.(theta.(phi))*dtheta
	dt(phi) = (den_t.(phi).*dnum_t.(phi) - 
            num_t.(phi).*dden_t.(phi))./
            (den_t.(phi).^2.0)
	dPsi_dt(t) = 1.0/pi./(1 .+ t.*t)
	dPsi(phi) = dPsi_dt(t(phi)).*dt(phi)
	x = reshape(u[:,1],1,m)
    dPsix = dPsi(x)
	dTu_u =  reshape([2.0 .+ dPsix;
    		  1.0 .+ dPsix; 
			  ones(m)';
			  ones(m)'], d, d, m)
    return dTu_u
end
function dCat_fd(u::Array{Float64,2}, s::Array{Float64,1}=[0.,0.])
	m, d = size(u)
	eps = 1.e-4

	u_p = copy(u)
	u_p[:,1] .+= eps
	u_m = copy(u)
	u_m[:,1] .-= eps
	dTx = mod.(Cat(u_p, s, 1)[2,:,:] - 
		   Cat(u_m, s, 1)[2,:,:],1)./(2*eps)
	
	u_p[:,1] .= u[:,1]
	u_p[:,2] .+= eps
	u_m[:,1] .= u[:,1] 
	u_m[:,2] .-= eps
	dTy = mod.(Cat(u_p, s, 1)[2,:,:] - 
		   Cat(u_m, s, 1)[2,:,:],1)./(2*eps)

	dTu_fd = reshape([dTx; dTy],d,d,m)
	return dTu_fd
end
