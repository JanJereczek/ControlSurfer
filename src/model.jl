#-----------------------------------------------------
# fixed CarbonCycle+Aerosol parameters
#-----------------------------------------------------
k_AU = 1/4              # 1/years
h_U = 150               # m. Thickness of upper ocean layer
h_D =  h_U*20           # m. Thickness of lower ocean layer
δ = h_D/h_U             # h_L/h_U dimensionless 
δDIC = 1.15             # dimensionless DIC_L/DIC_U at pre industrial
k_UD = (δ*δDIC)/1000    # 1/years
K0 = 0.0314806          # mol/(kg atm) Henry CO2 solubility
K1 = 1.32648E-6         # mol/kg
K2 = 9.19803E-10        # mol/kg
mA = 1.727E20           # moles in atmosphere
mO = 7.8E22             # moles in ocean
GtCtoppm(M_A) = (M_A * 1E6 * 1E15)/(12*mA) # ppm in atmosphere/ PgC in atmosphere
ppmtoGtC(conc) = conc*12*mA/(1E21)
W_U = mO*18E-3*h_U/(h_U+h_D) # whole upper ocean mass kg
W_D = mO*18E-3*h_D/(h_U+h_D) # whole lower ocean mass kg
pH_PI = 8.17
H_PI = 10^(-pH_PI) # mol/kg
CO2conc_a_PI = 280 # ppm
pCO2_a(CO2conc_a) = CO2conc_a*1E-6 # atm
Q = (K1/H_PI + 2*K1*K2/H_PI^2)*K0*pCO2_a(CO2conc_a_PI) # mol/kg
Qm = Q*W_U*12E-3*1E-12 # PgC
k_AL = 1/40            #1/yr
β_L = 1.7             # dimensionless [0.5 - 2.3]
#-----------------------------------------------------
# pre industrial initial conditions for the carbon reservoirs
#-----------------------------------------------------
M_A_PI = ppmtoGtC(CO2conc_a_PI) # PgC
M_U_PI = M_A_PI*(1 + K1/H_PI + K1*K2/H_PI^2)*W_U*K0/mA     # PgC
M_D_PI = M_U_PI*δ*δDIC    # PgC
M_L_PI = 2200 # PgC
#-----------------------------------------------------
# parameters entering temperature computation (greenhouse and heat transport parameters)
#-----------------------------------------------------
ECS = 3.5            # Equilibrium CliM_Ae Sensitivity. Celcius / doubling of CO2
TCR = 2.0            # Transient CliM_Ae Response. Celcius / doubling of CO2
F2X = 3.9            # Watts / m^2. Forcing due to a doubling of CO2 concentration
β = F2X/ECS          # Watts / (m^2 Celcius). Inverse equilibrium cliM_Ae sensitivity.
γ = F2X/TCR - β      # Watts / (m^2 Celcius). Thermal conductivity between layers
cvol = 0.13          # (Watts year) / (m^3 Celcius). Volumetric heat capacity of seawater.
TCR = F2X/(β+γ)      # Celsius per doubling of CO2 concentration
#-----------------------------------------------------
# parameters entering temperature computation (aerosol parameters taken from Helwegen2019, changed η)
#-----------------------------------------------------
#η = 0.742         # dimensionless
η = 1.0
αSO2 = 65        # Watts / m^2
βSO2 = 2246      # Mt of S / year
γSO2 = 0.23;     # dimensionless

#----------------------------------------------------
# a couple of auxiliary functions 
#----------------------------------------------------

function B(M_U) #dimensionles
    (sqrt( K1*( Qm*(K1-4K2)*(Qm-2M_U) + K1*M_U^2 ) ) - ( Qm*(K1-4K2) + M_U*(-K1+8K2) ) )/(2M_U*(K1-4K2))
end

function F(M_A,I) # antropogenic forcing
    if I > 0.0
        return F2X*log2(M_A/M_A_PI) - η*αSO2*exp(-(βSO2/I)^γSO2)
    else
        return F2X*log2(M_A/M_A_PI)
    end
end

function H(M_U) # Proton (hydrogen ion) concentration. mol/kg
    H = ( -Qm*K1 + K1*M_U + sqrt( K1*(Qm^2*(K1 - 4*K2) - 2*Qm*(K1 - 4*K2)*M_U + K1*M_U^2) ) )/(2*Qm)
    return H
end

function pH(M_U) # pH global scale
    return -log10(H(M_U))
end

function DIC(M_U) # μmol/kg
    return 1E6(M_U/12E-15)/(W_U)
end

function DIC_D(M_D) # μmol/kg
    return 1E6(M_D/12E-15)/(W_D)
end

function H2CO3(M_U) # μmol/kg
    H2CO3 = DIC(M_U)/(1 + K1/H(M_U) + K1*K2/H(M_U)^2)
    return H2CO3
end

function HCO3(M_U) # μmol/kg
    HCO3 = K1*H2CO3(M_U)/H(M_U)
    return HCO3
end

function CO3(M_U) # μmol/kg
    CO3 = K2*HCO3(M_U)/H(M_U)
    return CO3
end

CO3sat = CO3(M_U_PI)/3.44 

function Ω(M_U)
    return CO3(M_U)/CO3sat
end


#-----------------------------------------------------
# ICE parameters and functions
#-----------------------------------------------------

function Vmcons(model_parameters)
    
    Tp, Tm, Vp, Vm, τmelt, τfreeze = model_parameters
    
    x = ((-Tm + Tp)/(Tm + Tp + 2*sqrt(Tm*Tp)))^(1/3)
    Vm = ( -2 + Vp*(1 + x + 1/x ) )/( -1 + x + 1/x )
    
    return Vm
end

Greenland_params = [1.52, 0.3, 0.77, 0.3526554620064224, 470.0, 5000.0]#[Tp, Tm, Vp, Vm, tau_melt, tau_freeze] 
Greenland_params[4] = Vmcons(Greenland_params)
Antarctica_params = [6.8, 4, 0.44, 0.07857839308355193, 3000.0, 5500.0]#[Tp, Tm, Vp, Vm, tau_melt, tau_freeze] 
Antarctica_params[4] = Vmcons(Antarctica_params)
# sea level rise potential
SLRpotentialG = 7.4
SLRpotentialA = 55



function dV_dt(V, Tf, model_parameters)

    Tp, Tm, Vp, Vm, τmelt, τfreeze = model_parameters
    
    a = 3*(Vm + Vp)/2
    b = -3*Vm*Vp
    c = (Vp - Vm)^3/(2*(Tm - Tp))
    d = ( Tp*Vm^2*(Vm-3Vp) - Tm*Vp^2*(Vp-3Vm) )/(2*(Tm - Tp))
    
    function μ(V,Tf)
        if (- V^3 + a*V^2 + b*V + c*Tf + d) > 0
            return 1/τfreeze
        else
            if V < 1.0e-4
                return 0.0
            else
                return 1/τmelt
            end
        end
    end
       
    return μ(V,Tf)*(- V^3 + a*V^2 + b*V + c*Tf + d)
end

dV_dtA(V,Tf) = dV_dt(V,Tf,Antarctica_params)
dV_dtG(V,Tf) = dV_dt(V,Tf,Greenland_params) 


#-----------------------------------------------------------------------------------
# Sea level rise
#-----------------------------------------------------------------------------------
# SLR ice
SLR_G(VG) = SLRpotentialG*(1-VG)
SLR_A(VA) = SLRpotentialA*(1-VA)
SLRice(VG,VA) =  SLR_G(VG) + SLR_A(VA)
rateSLRice(VG,VA,δT_U) = - SLRpotentialG*(dV_dtG(VG,δT_U))*100 - SLRpotentialA*(dV_dtA(VA,δT_U))*100

# SLR thermal
αU = 2.3e-4 #K^-1
αD = 1.3e-4 #K^-1
SLRthermal(δT_U,δT_D;αU = αU, αD = αD) = αU*δT_U*h_U + αD*δT_D*h_D

# SLR glacier
τ_gl = 200 # years
SLRpotential_gl = 0.5 # meters
ζ = 2 # Related to equlibrium SLRglacier sensitivity to temperature change (Celsius)
SLRgl_eq(δT_U) = SLRpotential_gl*tanh(δT_U/ζ)
SLRglacier(S) = S

# SLRtotal(δT_U,δT_D,VG,VA,S)
SLRtotal(δT_U,δT_D,VG,VA,S) = SLRthermal(δT_U,δT_D) + SLRice(VG,VA) + SLRglacier(S)

function model!(du,u,p,t)
    M_A, M_U, M_D, M_L, δT_U, δT_D, VG, VA, S = u
    
    Emissions, Injections = p
    
    du[1] = Emissions(t) - k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) - k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))
    du[2] = k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) - k_UD*(M_U - M_D/(δ*δDIC))
    du[3] = k_UD*(M_U - M_D/(δ*δDIC))
    du[4] = k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))
#     du[5] = dδT_U = ( F(M_A,Injections(t)) - β*δT_U - γ*(δT_U - δT_D) )/(cvol*h_U)
    du[5] = ( F2X*log2(M_A/M_A_PI) - β*δT_U - γ*(δT_U - δT_D) )/(cvol*h_U)
    du[6] = γ*(δT_U - δT_D)/(cvol*h_D)
    du[7] = dV_dtG(VG,δT_U)
    du[8] = dV_dtA(VA,δT_U)
    du[9] = (1/τ_gl)*(SLRgl_eq(δT_U) - S)
     
end

function controlled_model(u,p,t)
    du = similar(u)
    uw, K = p["yw"], p["K"]
    M_A, M_U, M_D, M_L, δT_U, δT_D, VG, VA, S = u[1:9]
    x_control = u[10:end]
    e = uw - u[1]
    control = K.C * x_control + K.D .* e

    du[1] = control[1] - k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) - k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))
    du[2] = k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) - k_UD*(M_U - M_D/(δ*δDIC))
    du[3] = k_UD*(M_U - M_D/(δ*δDIC))
    du[4] = k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))
    du[5] = ( F2X*log2(M_A/M_A_PI) - β*δT_U - γ*(δT_U - δT_D) )/(cvol*h_U)
    du[6] = γ*(δT_U - δT_D)/(cvol*h_D)
    du[7] = dV_dtG(VG,δT_U)
    du[8] = dV_dtA(VA,δT_U)
    du[9] = (1/τ_gl)*(SLRgl_eq(δT_U) - S)
    du[10:end] = vec(K.A * x_control + K.B .* e)
    return du
end