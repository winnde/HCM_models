export Plots, LinearAlgebra
using Printf

#=
This code implements Hydro-Mechanical model II : Nonlinear diffusion equation for Porosity from Schmalholz et al. 2024
Equation of porosity diffusion with porosity-dependent permeability 
Without (de)hydratation reactions
Including poroelastic effects

Equations used to find the PDE of porosity in this case : 
1- Fluid mass conservation (per unit volume)
  ∂/∂t(ϕρf) = ∂/∂x(vfϕρf)       [1]
2- Dracy's law
  ϕvf = k/ηf ∂Pf/dx             [2]
3- Pore compressibility equation (poroelasticity, Yarushina & Podladchikov, 2015)
  ∂ϕ/∂t = βϕ ∂Pf/∂t             [3]
4 - Permeability–porosity relationship (Kozeny–Carman formulation)
  k = kref (ϕ/ϕref)ⁿ            [7]

Resulting nonlinear porosity diffusion equation:
∂ϕ/∂t = ∂/∂x (kref/(Βϕ*ηf) (ϕ/ϕref)^n ∂ϕ/dx)   [8]
=#

function implicit_pt_HM_II(;nt,nvx, tol=1e-6, niter=1e5)
 
  # Part A : Definition of parameters
  # ----------------------------------

  # Physical parameters
  L_phys          = 1          # Length                             | [m]
  ηf_phys         = 1e-3       # Fluid viscosity                    | [Pa.s]
  Βϕ_phys         = 0.1        # Pore compressibility               | [Pa-1]
  kmin_phys       = 0.1        # Minimum initial permeability       | [m2] 
  kmax_phys       = 1          # Maximum initial permeability       | [m2]
  kref_phys       = 0.2        # Reference permeability             | [m2]
  ϕ_amb           = 0.1        # Ambiant porosity 
  dϕ_phys         = 0.1        # Magnitude of ϕ perturbation        | [m] -> For initial condition
  ϕ_max           = 0.8        # Maximum of porosity anomaly              -> For initial condition
  n_darcy         = 3          # Carman Kozeny exponent
  ϕref            = 0.2

  # Characteristic parameters
  Lc          = L_phys                 # Lengh scale [m]
  tc          = ηf_phys*Βϕ_phys        # time scale  [s]
  ηc          = ηf_phys                # viscosity scale [Pa.s]

  # Nondimensional properties
  Lx          = L_phys  / Lc
  ηf          = ηf_phys / ηc
  Βϕ          = Βϕ_phys / (ηc/tc)
  kmin        = kmin_phys  / Lc^2
  kmax        = kmax_phys  / Lc^2
  kref        = kref_phys / Lc^2 
  dϕ          = dϕ_phys / Lc
  
  κ           = kref/Βϕ/ηf              # Coeff. diff.

  # Numerics
  ncx         = nvx-1                   # Numerical resolution in centroids
  Δx          = Lx/ncx
  xv          = LinRange(-Lx/2, Lx/2, nvx)
  xc          = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, ncx)
  
  Δt          = Δx^2/(2.1*κ) 
  θ           = 4.3* (Lx)/nvx 

  # Allocate arrays 
  ϕ0          = ϕ_amb * ones(nvx)
  ϕ           = zero(ϕ0)
  ϕ_init      = zero(ϕ0)
  k0          = zero(ϕ0)
  qϕ          = zeros(ncx)
  div_qϕ      = zeros(ncx-1)
  R_ϕ         = zeros(ncx-1)
  ∂ϕ∂τ        = zeros(ncx-1)
  coeff_diff  = zero(qϕ)


  # Part B : Initialization
  # -----------------------

  # Initial condition for porosity :
  ϕ0              .= ϕ0 + ϕ_max.*exp.(-xv.^2/2/dϕ^2)
  ϕ                = copy(ϕ0)
  ϕ_init           = copy(ϕ0)

  # Initial condition for permeability
  k0               = kmin * ones(length(k0))
  for i=1:ncx
    xv[i] > -1/5*Lx && xv[i] < 1/5*Lx ? k0[i] = kmax : nothing
  end


  # Part C : Solve equation
  # -----------------------

  # Time loop
  for it=1:nt 
      
    @printf("Time step = %06d\n", it)
      
    #update
    ϕ0    .= ϕ
    # Pseudo-time step calculation
    Δτ     = 1e6 * Δt * Δx^2 / ( maximum(kref/(Βϕ*ηf) * (ϕ./ϕref).^n_darcy) * 4)      # Don't ask me why we use 1e6...🥲

    # Pseudo-time loop
    for iter=1:niter

      # Limits (v=0)
      ϕ[1]          =  ϕ[2]
      ϕ[end]        =  ϕ[end-1]              
      
      # Diffusion factor calculation
      # We make interpolation to have a good indice number and to be in phase with qϕ
      coeff_diff    = @. (kref/(Βϕ*ηf) * (ϕ[1:end-1]/ϕref)^n_darcy + kref/(Βϕ*ηf) * (ϕ[2:end]/ϕref)^n_darcy )/2
      
      # flux
      qϕ           .= - coeff_diff .* diff(ϕ) /Δx
      
      # ∇q
      div_qϕ       .= diff(qϕ) ./ Δx

      # Residuals 
      R_ϕ          .= - ((ϕ[2:end-1]  .- ϕ0[2:end-1])  ./ Δt + div_qϕ)
      
      # Pseudo rate update
      ∂ϕ∂τ         .= R_ϕ .+ (1.0-θ).*∂ϕ∂τ    

      # Solution update 
      ϕ[2:end-1]  .+= Δτ  .* ∂ϕ∂τ 

      # Error calculation
      errϕ          = norm(R_ϕ) / length(R_ϕ)
      if mod(iter, 1000)==0
          @printf("iter = %06d:  Err. ϕ = %1.6e\n", iter, errϕ)
      end 
      errϕ<tol ? break : nothing
      iter == niter ? error("That don't converge !") : nothing

    end 

    # Part D : Visualisation
    # ----------------------

    p1 = plot(xv, ϕ0, label="ϕ0", title="Porosity-dependent permeability, step $it ", xlabel="x", ylabel="ϕ")
    p1 = plot!(xv, k0, label="k_init")
    p1 = plot!(xv, ϕ_init, label="ϕ_init")
    p1 = plot!(xv, ϕ, label="ϕ")
    display(plot(p1))
    sleep(0.05)
      
  end
 
end

implicit_pt_HM_II(;nt=500, nvx=151, tol=1e-6, niter=1e5)