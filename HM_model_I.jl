using Plots, SpecialFunctions, LinearAlgebra
using Printf

# Explanation:
# -------------

# This code solves the PDEs of Model I from Schmalholz et al., 2024:
# Equations for porosity and Pf (fluid pressure) diffusion
# without (de)hydration reactions,
# and including poroelastic effects.

# PDEs:
# ∂Pf/∂t = ∂/∂x (k / (B * ϕ * ηf) * ∂Pf/∂x)   [1] Eq. 5 in Schmalholz et al., 2024
# ∂ϕ/∂t  = ∂/∂x (k / (B * ϕ * ηf) * ∂ϕ/∂x)    [2] Eq. 6 in Schmalholz et al., 2024

# Four numerical methods are proposed to solve these equations:
# explicit_HM_I()               = explicit method
# pt_method_no_damp()           = implicit pseudo-transient method without damping
# pt_method()                   = implicit pseudo-transient method with damping
# classical_implicit_method()   = implicit method (matrix inversion)

function main(;nt, nvx)

    # Which method do you want?
    # -------------------------
    explicit            = 0
    pt_no_damp          = 0
    pt                  = 1
    classical_implicit  = 0

    Methods = explicit + pt_no_damp + pt + classical_implicit
    Methods > 1 ? error("You can’t choose more than one method!") : nothing
    Methods == 0 ? error("Please, select a method...") : nothing


    # I - Initialization
    #--------------------

    # Physical parameters 
    L_phys      = 1             # Length of the box             | [m]
    ηf_phys     = 1000          # Viscosity of fluid            | [Pa.s]
    βϕ_phys     = 0.1           # Pore compressibility          | [Pa-1]
    k_phys      = 0.5           # Permeability                  | [m2]
    P_BG_phys   = 1.6e9         # Pressure background           | [Pa]
    dPf_phys    = 0.1           # Magnitude of Pf pertubation   | [m]
    Pf_max_phys = 0.8e9         # Maximum of Pf anomaly         | [Pa]
    ϕ_amb       = 0.1           # Ambiant porosity 
    dϕ_phys     = 0.1           # Magnitude of ϕ perturbation   | [m]
    ϕ_max       = 0.5           # Maximum of porosity anomaly
    

    # Characteristic parameters
    Lc          = L_phys                # LenghScale      | [m]
    Pc          = P_BG_phys             # Pscale          | [Pa]
    tc          = ηf_phys/βϕ_phys       # timescale       | [s]
    ηc          = ηf_phys               # viscosity scale | [Pa.s]

    # Nondimensional properties
    Lx          = L_phys / Lc 
    P_BG        = P_BG_phys / Pc
    ηf          = ηf_phys / ηc
    Βϕ          = βϕ_phys * Pc
    k           = k_phys  / Lc^2 

    dϕ          = dϕ_phys / Lc       # for initial conditions
    dPf         = dPf_phys / Lc 
    Pf_max      = Pf_max_phys / Pc

    # Mesh setup and Numeric
    ncx         = nvx-1        # numerical resolution in centroids
    Δx          = Lx/ncx
    xv          = LinRange(-Lx/2, Lx/2, nvx)
    xc          = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, ncx)
  
    # Time and pseudo-time
    κ           = k/Βϕ/ηf                           # Diffusivity coeff
    Δt          = Δx^2/(2.1*κ)                      # For stabily of explicit method, Δt need to be <1/2 (it's why there are 2.1)
    Β           = κ * Δt / Δx^2                     # For classical implicit method 
    Δτ          = 2 * Δt*Δx^2/(4.1*κ*Δt + Δx^2)     # For pseudo-transient method
    θ           = Lx/nvx * 4.3                      # For damping

    # Allocate arrays
    ϕ           = zeros(nvx)
    ϕ0          = zero(ϕ)
    ϕ_init      = zero(ϕ)
    qϕ          = zeros(ncx)
    div_qϕ      = zeros(ncx-1)
    R_ϕ         = zeros(ncx-1)
    ∂ϕ∂τ        = zeros(ncx-1)

    Pf0         = zeros(nvx)
    Pf          = zeros(nvx)
    Pf_init     = zeros(nvx)
    qPf         = zeros(ncx)
    div_qPf     = zeros(ncx-1)
    R_Pf        = zeros(ncx-1)
    ∂Pf∂τ       = zeros(ncx-1)

    if classical_implicit == 1
        #if @isdefined classical_implicit_method
        
        # Definition of matrix : 
        Pf_centre        = ones(nvx) * (1 - 2*Β)     # Creation de la diagonale centre 
        Pf_centre[1]     = 1                         # Condition aux limites aux bord inférieur -> bords fixes
        Pf_centre[end]   = 1                         # Condition aux limites aux bord superieur -> Bord fixes
        Pf_sup           = ones(nvx-1) * Β           # Création de la diagonal supérieur 
        Pf_sup[1]        = 0
        Pf_inf           = ones(nvx-1) * Β           # Creation de la diagonal inférieur 
        Pf_inf[end]      = 0
        mat_Pf           = diagm(0 => Pf_centre, 1 => Pf_sup, -1 => Pf_inf) # creation de la diagonal tridiagonalisée

        ϕ_centre         = ones(nvx) * (1 - 2*Β)     # Creation de la diagonale centre 
        ϕ_centre[1]      = 1                         # Condition aux limites aux bord inférieur 
        ϕ_centre[end]    = 1                         # Condition aux limites aux bord superieur 
        ϕ_sup            = ones(nvx-1) * Β           # Création de la diagonal supérieur 
        ϕ_sup[1]         = 0
        ϕ_inf            = ones(nvx-1) * Β           # Creation de la diagonal inférieur 
        ϕ_inf[end]       = 0
        mat_ϕ            = diagm(0 => Pf_centre, 1 => Pf_sup, -1 => Pf_inf) # creation de la diagonal tridiagonalisée

    end


    # II - Initialization
    # --------------------
  
    # Initial condition for porosity :
    ϕ                = ϕ_amb * ones(nvx)
    ϕ               .= ϕ + ϕ_max.*exp.(-xv.^2/2/dϕ^2)
    ϕ_init           = copy(ϕ)

    # Initial condition for Pf :
    Pf               = P_BG * ones(nvx)
    Pf              .= Pf + Pf_max.*exp.(-xv.^2/2/dPf^2)
    Pf_init          = copy(Pf)

        
    # III - Computation
    # ------------------

    # Time loop
    for it=1:nt 
        
        @printf("Time step = %06d\n", it)

        # update 
        Pf0     .= Pf 
        ϕ0      .= ϕ

        if explicit == 1
            explicit_method(;Pf, ϕ, Pf0, ϕ0, qϕ, qPf, κ, Δx, Δt)
        end 

        if pt_no_damp == 1 
            pt_method_no_damp(;tol=1e-6, niter=1e5, Pf, ϕ, qϕ, qPf, Δx, Δt, div_qPf, div_qϕ, 
            R_Pf, R_ϕ, Pf0, ϕ0, k, Βϕ, ηf, Δτ)
        end 
        
        if pt == 1 
            pt_method(;tol=1e-6, niter=1e5, Pf, ϕ, qϕ, qPf, κ, Δx, Δt, div_qPf, div_qϕ, 
            R_Pf, R_ϕ, ∂Pf∂τ, ∂ϕ∂τ, k, Βϕ, ηf, Pf0, ϕ0, Δτ, θ)
        end 

        if classical_implicit == 1 
            classical_implicit_method(;Pf, mat_Pf, ϕ, mat_ϕ)
        end


        # VI - Visualization 
        # -------------------

        # Analytical solution
        t            = it*Δt
        ϕana         = zeros(nvx)
        Pfana        = zeros(nvx)
        ϕana        .= @. ϕ_amb + ϕ_max  ./ sqrt.(1 + 2*t*κ/dϕ^2) .* exp.(-xv.^2. / (2*(dϕ^2 + 2*t*κ)) )
        Pfana       .= @. P_BG + Pf_max ./ sqrt.(1 + 2*t*κ/dPf^2) .* exp.(-xv.^2. / (2*(dPf^2 + 2*t*κ)) )

        # Plots
        p1 = plot(xv, Pf_init,label="Pf init", title = "Pf", xlabel="x", ylabel="Pf")
        p1 = scatter!(xv, Pf0, label="Pf0", markercolor=:white, markerstrokecolor=:black)
        p1 = scatter!(xv, Pf, label="Pf",  markercolor=:white, markerstrokecolor=:blue)
        p1 = plot!(xv, Pfana, label="Pf ana")
        p2 = plot(xv, ϕ_init, label="ϕ init", title = "ϕ", xlabel="x", ylabel="Pf")
        p2 = scatter!(xv, ϕ0, label="ϕ0", markercolor=:white, markerstrokecolor=:black)
        p2 = scatter!(xv, ϕ, label="ϕ", markercolor=:white, markerstrokecolor=:blue)
        p2 = plot!(xv, ϕana, label="ϕ ana")
        
        display(plot(p1,p2))
        sleep(0.1)
        
    end

end


function explicit_method(;Pf, ϕ, Pf0, ϕ0, qϕ, qPf, κ, Δx, Δt)

    # limits -> v = 0
    Pf0[1]          =  Pf0[2]
    Pf0[end]        =  Pf0[end-1]
    ϕ0[1]           =  ϕ0[2]
    ϕ0[end]         =  ϕ0[end-1]

    #flux : 
    qϕ             .= - κ.* diff(ϕ0)/Δx
    qPf            .= - κ.* diff(Pf0)/Δx

    # Explicit update
    Pf[2:end-1]   .-= diff(qPf)/Δx .* Δt
    ϕ[2:end-1]    .-= diff(qϕ) /Δx .* Δt
    
end

function pt_method_no_damp(;tol=1e-6, niter=1e5, Pf, ϕ, qϕ, qPf, Δx, Δt, div_qPf, div_qϕ, 
    R_Pf, R_ϕ, Pf0, ϕ0, k, Βϕ, ηf, Δτ)

    # Pseudo-time loop
    for iter=1:niter
      
        # limits -> v=0
        Pf[1]            =  Pf[2]
        Pf[end]          =  Pf[end-1]
        ϕ[1]             =  ϕ[2]
        ϕ[end]           =  ϕ[end-1]
        
        #flux 
        qϕ              .= - k .* diff(ϕ)/Δx
        qPf             .= - k .* diff(Pf) ./Δx
  
        # ∇q 
        div_qPf         .=  diff(qPf) ./ Δx
        div_qϕ          .=  diff(qϕ)  ./ Δx
  
        # Residuals 
        R_Pf            .= - ((Βϕ*ηf)*(Pf[2:end-1] .- Pf0[2:end-1]) ./ Δt + div_qPf )
        R_ϕ             .= - ((Βϕ*ηf)*(ϕ[2:end-1]  .- ϕ0[2:end-1])  ./ Δt + div_qϕ)
        
        # Solution update
        Pf[2:end-1]  .+= Δτ/(ηf*Βϕ)  .* R_Pf 
        ϕ[2:end-1]   .+= Δτ/(ηf*Βϕ)  .* R_ϕ   
  
        # Error calculation
        errPf         = norm(R_Pf) / length(R_Pf)
        errϕ          = norm(R_ϕ)  / length(R_ϕ) 
        if mod(iter, 10)==0
            @printf("iter = %06d:  Err. Pf = %2.1e --- Err. ϕ = %2.1e\n", iter, errPf, errϕ)
        end
        errPf<tol ? break : nothing
  
    end 


end

function pt_method(;tol=1e-6, niter=1e5, Pf, ϕ, qϕ, qPf, κ, Δx, Δt, div_qPf, div_qϕ, 
    R_Pf, R_ϕ, ∂Pf∂τ, ∂ϕ∂τ, k, Βϕ, ηf, Pf0, ϕ0, Δτ, θ)

    for iter=1:niter
      
        # Limits -> v=0
        Pf[1]            =  Pf[2]
        Pf[end]          =  Pf[end-1]
        ϕ[1]             =  ϕ[2]
        ϕ[end]           =  ϕ[end-1]
        
        #flux
        qPf             .= - k .* diff(Pf) ./Δx
        qϕ              .= - k .* diff(ϕ)/Δx
  
        # ∇q
        div_qPf         .=  diff(qPf) ./ Δx
        div_qϕ          .=  diff(qϕ)  ./ Δx
  
        # Residuals 
        R_Pf            .= - ((Βϕ*ηf)*(Pf[2:end-1] .- Pf0[2:end-1]) ./ Δt + div_qPf )
        R_ϕ             .= - ((Βϕ*ηf)*(ϕ[2:end-1]  .- ϕ0[2:end-1])  ./ Δt + div_qϕ)
        
        # Pseudo rate update
        ∂Pf∂τ        .= R_Pf .+ (1.0-θ) .* ∂Pf∂τ 
        ∂ϕ∂τ         .= R_ϕ .+ (1.0-θ)  .* ∂ϕ∂τ 

        # Solution update
        Pf[2:end-1]  .+= Δτ/(ηf*Βϕ)  .* ∂Pf∂τ 
        ϕ[2:end-1]   .+= Δτ/(ηf*Βϕ)  .* ∂ϕ∂τ   
  
        # Error calculation
        errPf         = norm(R_Pf) / length(R_Pf)
        errϕ          = norm(R_ϕ)  / length(R_ϕ)
        if mod(iter, 100)==0
            @printf("iter = %06d:  Err. Pf = %2.1e --- Err. ϕ = %2.1e\n", iter, errPf, errϕ)
        end
        errPf<tol ? break : nothing
  
    end 

end

function classical_implicit_method(;Pf, mat_Pf, ϕ, mat_ϕ)
    
    Pf  .= transpose.(mat_Pf * Pf)
    ϕ   .= transpose.(mat_ϕ * ϕ)

end


main(;nt=10, nvx=101)