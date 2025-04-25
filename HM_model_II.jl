export Plots, LinearAlgebra
using Printf

#=
This code makes Hydro-Mechanical Model II : Nonlinear diffusion equation for Porosity from Schmalholz et al. 2024
  - Équation of porosity with porosity-dependent permeability 
  - Without reaction of (de)hydratation
  - With a poroelasticity 

Equations used to find the PDE of porosity in this case : 
1- Eq for the conservation of fluid mass per unit volume
  ∂/∂t(ϕρf) = ∂/∂x(vfϕρf)       [1]
2- Dracy law
  ϕvf = k/ηf ∂Pf/dx             [2]
3- Eq of pore compressibility in classical poroelasticity (Yarushina & Poldlachikov, 2015)
  ∂ϕ/∂t = βϕ ∂Pf/∂t             [3]
4 - Eq that give the nonlinar permeability-porosity relationship (Kozeny Carman formulation)
  k = kref (ϕ/ϕref)ⁿ            [7]

Give the PDE of diffusion of porosity : 
∂ϕ/∂t = ∂/∂x (kref/(Βϕ*ηf) (ϕ/ϕref)^n ∂ϕ/dx)   [8]
=#

function implicit_pt_HM_II(;nt,nvx, tol=1e-6, niter=1e5)
  
 
  # Part A : Definition of parameters
  # ----------------------------------

  # Physical parameters
  L_phys          = 1          # Dimention                          | [m]
  ηf_phys         = 1e-3       # Fluid viscosity                    | [Pa.s]
  Βϕ_phys         = 0.1        # Compressibilité des pores          | [Pa-1]
  kmin_phys       = 0.1        # minimum initial permeability       | [m2] 
  kmax_phys       = 1          # maximum intial permeability        | [m2]
  kref_phys       = 0.2        # reference permeability             | [m2]
  ϕ_amb           = 0.1        # ambiant porosity 
  dϕ_phys         = 0.1        # Magnitude of ϕ perturbation        | [m] -> For initial condition
  ϕ_max           = 0.8        # Maximum of porosity anomaly              -> For initial condition
  n_darcy         = 3          # Carman Kozeny exponent
  ϕref            = 0.2

  # Characteristic parameters
  Lc          = L_phys                 # LenghScale [m]
  tc          = ηf_phys*Βϕ_phys        # timescale  [s]
  ηc          = ηf_phys                # viscosity scale [Pa.s]

  # Nondimensional properties
  Lx          = L_phys  / Lc
  ηf          = ηf_phys / ηc
  Βϕ          = Βϕ_phys / (ηc/tc)
  kmin        = kmin_phys  / Lc^2
  kmax        = kmax_phys  / Lc^2
  kref        = kref_phys / Lc^2 
  dϕ          = dϕ_phys / Lc
  κ           = kref/Βϕ/ηf

  # Numerics
  ncx         = nvx-1               # numerical resolution in centroids
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

  #========================
  Part B : Initialization
  ========================#

  # Initial condition for porosity :
  ϕ0              .= ϕ0 + ϕ_max.*exp.(-xv.^2/2/dϕ^2)
  ϕ                = copy(ϕ0)
  ϕ_init           = copy(ϕ0)

  # Initial condition of permeability
  k0               = kmin * ones(length(k0))
  
  for i=1:ncx
    xv[i] > -1/5*Lx && xv[i] < 1/5*Lx ? k0[i] = kmax : nothing
  end

  #=========================
  Part C : Computation
  =========================#

  for it=1:nt 
      
    @printf("Time step = %06d\n", it)
      
    #update
    ϕ0    .= ϕ
    Δτ     = 1e6 * Δt * Δx^2 / ( maximum(kref/(Βϕ*ηf) * (ϕ./ϕref).^n_darcy) * 4)

    # Pseudo-time loop (iterations)
    for iter=1:niter

      # Definition des bords (V = 0)
      ϕ[1]        =  ϕ[2]
      ϕ[end]      =  ϕ[end-1]              
      
      # diffusion factor calculation : 
      coeff_diff     .= @. (kref/(Βϕ*ηf) * (ϕ[1:end-1]./ϕref).^n_darcy + kref/(Βϕ*ηf) * (ϕ[2:end]./ϕref).^n_darcy )/2
      # On fait une moyenne de ϕ pour avoir le bon nombre d'indice pour coeff_diff
      
      # flux
      qϕ             .= - coeff_diff .* diff(ϕ) /Δx
      
      # ∇q
      div_qϕ         .= diff(qϕ) ./ Δx

      # Residuals 
      R_ϕ            .= - ((ϕ[2:end-1]  .- ϕ0[2:end-1])  ./ Δt + div_qϕ)
      
      # Pseudo rate update
      ∂ϕ∂τ           .= R_ϕ .+ (1.0-θ).*∂ϕ∂τ    

      # Update 
      ϕ[2:end-1]    .+= Δτ  .* ∂ϕ∂τ 

      # Calculation of error
      errϕ            = norm(R_ϕ) / length(R_ϕ)
      if mod(iter, 100)==0
          @printf("iter = %06d:  Err. ϕ = %1.6e\n", iter, errϕ)
      end 
      errϕ<tol ? break : nothing
      iter == niter ? error("That don't converge !") : nothing

    end 

    #=========================
      Part D : Visualisation
    =========================#

    p1 = plot(xv, ϕ0, label="ϕ0", title="Porosity-dependent permeability, $it ")
    p1 = plot!(xv, k0, label="k_init")
    p1 = plot!(xv, ϕ_init, label="ϕ_init")
    p1 = plot!(xv, ϕ, label="ϕ")
    display(plot(p1))
    sleep(0.05)
      
  end
 
end

implicit_pt_HM_II(;nt=500, nvx=151, tol=1e-6, niter=1e5)