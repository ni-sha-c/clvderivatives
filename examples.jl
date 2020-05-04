using LinearAlgebra
"""
	n: number of timesteps
    m: number of initial conditions
    s[0] = abs(lambda), s[1] = alpha
    output size: mxdx(n+1)
"""
function step(u::Array{Float64,2}, s::Array{Float64,1}
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
function inverse_step(u::Array{Float64,2}, s::Array{Float64,1}=[0.,0.], n::Int64=1)
	
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
    return u_trj
end
"""
	Compute Lyapunov Exponents and Covariant 
	Lyapunov vectors.

"""
function clvs(DTu::Array{Float64,3},du::Int64)
	d = size(DTu)[1]
	m = size(DTu)[3]
	lyap_exps = zeros(du)
	R = zeros(du,du,m)
	Q = zeros(d,du,m)
	Q[:,:,1], R[:,:,1] = qr!(rand(d,du))
	for i=2:m
		Q[:,:,i] = DTu[:,:,i-1]*Q[:,:,i-1]
		Q[:,:,i], R[:,:,i] = qr!(Q[:,:,i])
		lyap_exps .+= log.(abs.(diag(R[:,:,i])))./m
	end
	C = zeros(du,du,m)
	C[:,:,end] = diagm(ones(du))	
	for i=reverse(1:m-1)
		C[:,:,i] = R[:,:,i+1]\C[:,:,i+1]
		[normalize!(view(C,:,j,i)) for j = 1:du]
		Q[:,:,i] = Q[:,:,i]*C[:,:,i]
	end
	return lyap_exps, Q
end
function dstep(u::Array{Float64,2}, s::Array{Float64,1}=[0.,0.])
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

    dPsix = dPsi(u[:,1])'
	dTu_u =  reshape([2.0 .+ dPsix;
			  ones(m)';
    		  1.0 .+ dPsix; 
			  ones(m)'], d, d, m)
    return dTu_u
end

u = rand(1,2)
s = [0.7,0.3]
n = 1000
un = step(u, s, n)[:,:,1]
dTu = dstep(un, s)
les, covlyapvecs = clvs(dTu,2)

