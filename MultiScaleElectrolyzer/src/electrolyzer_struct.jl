struct Electrolyzer
    # Device Parameters
    ϕ::Float64 # Nominal capacity of electrolyzer
    α_max::Float64 # maximum production coef of electrolyzer (kg/MWh)
    α_min::Float64 # minimum production coef of electroylzer (kg/MWh)
    ξ_start::Float64 # electrical efficiency of new electrolyzer (% LHV)
    η_overpotential::Float64 # overpotential [μV/hr]
    V_i::Float64 # voltage for each cell [V]
    ℓ::Float64 # lifetime of the electrolyzer
    δ_on::Float64 # Normal degradation coef
    η_startup::Float64 # Start degradation overpotential [μV/hr]
    δ_start::Float64 # startup degradation coef
    ρ_sb::Float64

    # Economics
    λ_H::Float64 # cost of hydrogen $/kg
    λ_OPEX::Float64 # cost of operation of stack $/kg
    λ_CAPEX_Plant::Float64 # cost of stack $
    λ_CAPEX_Stack::Float64 # cost of stack $
    i::Float64 # discount rate %
    CRF::Float64 # capital recovery factor $/MW
    λ_CAPEX_yearly::Float64 # yearly CAPEX cost
    λ_CAPEX_replace::Float64 # cost to replace a stack


    function Electrolyzer(ϕ, ξ_start, η_overpotential, V_i, ℓ, η_startup, λ_H, i, λ_CAPEX_Plant, λ_CAPEX_Stack, ρ_sb)
        LHV = 33.36 # lower heating value (or net heating value) for hydrogen

        HHV = 39.4 # kWh/kg, higher heating value for hydrogen
        V_o = 1.48 # thermoneutral voltage


        α_max = 1000*ξ_start/LHV # calculating the maximum production coef for the system

        # δ_on = (η_overpotential/V_i)/(1000*LHV) # calculating the decrease in production per hour of use

        δ_on = (V_o * 1000)/HHV * (η_overpotential/(V_i * (V_i + η_overpotential)))

        α_min = α_max - δ_on * ℓ * 8760


        # δ_start = 1000*(η_startup)/HHV
        # δ_start = (V_o * 1000)/HHV * ((0.0204/500)/(V_i * (V_i + (0.0204/500))))
        δ_start = (V_o * 1000)/HHV * ((η_startup)/(V_i * (V_i + (η_startup))))
        # derived from voltage increase per 500 reactions running at 0.6 A/cm^2 per the DOE targets

        # δ_start = ((0.0204/500)/V_i)/(1000*LHV)


        # λ_CAPEX = ϕ * 248 * 1000

        CRF = (i * (1 + i)^ℓ)/((1 + i)^ℓ - 1)
        λ_CAPEX_yearly = CRF * λ_CAPEX_Plant * ϕ
        λ_OPEX = λ_CAPEX_yearly * 0.02 # H. Nami et al 2022

        λ_CAPEX_replace =  λ_CAPEX_Stack * ϕ # H. Nami et al 2022 TODO: Check this figure and understand it better

        new(ϕ, α_max, α_min, ξ_start, η_overpotential, V_i, ℓ, δ_on, η_startup, δ_start, ρ_sb, λ_H, λ_OPEX, λ_CAPEX_Plant, λ_CAPEX_Stack, i, CRF, λ_CAPEX_yearly, λ_CAPEX_replace)
    end
end
Base.show(io::IO, x::Electrolyzer) = print(io, """
Electrolyzer Struct:
  Nominal Capacity (ϕ): $(x.ϕ) MW
  Production Coefficient: α_min = $(x.α_min) kg/MWh, α_max = $(x.α_max) kg/MWh
  Efficiency and Degradation:
    Initial Efficiency (ξ_start): $(x.ξ_start*100) % LHV
    Overpotential Rate (η): $(x.η_overpotential) μV/hr
    Start-up efficiency loss (η): $(x.η_startup) %
    Voltage per Cell (Vᵢ): $(x.V_i) V
    Lifetime (ℓ): $(x.ℓ) yrs
    Degradation Coeffs: δ_on = $(x.δ_on), δ_start = $(x.δ_start)
  Economics:
    Hydrogen Price (λ_H): \$$(x.λ_H)/kg
    OPEX (λ_OPEX): \$$(x.λ_OPEX) /yr
    CAPEX (λ_CAPEX_Plant): \$$(x.λ_CAPEX_Plant) /MW
    CAPEX (λ_CAPEX_Stack): \$$(x.λ_CAPEX_Stack) /MW
    Discount Rate (i): $(x.i*100) %
    CRF: $(x.CRF)
    Yearly CAPEX (λ_CAPEX_yearly): \$$(x.λ_CAPEX_yearly)/year
    Replacement Cost (λ_CAPEX_replace): \$$(x.λ_CAPEX_replace)
""")
