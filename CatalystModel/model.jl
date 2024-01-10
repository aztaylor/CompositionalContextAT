using Catalyst, ModelingToolkit, Plots, Latexify, Main.UnitfulBio, Unitful, IfElse, UnitfulAngles
# Have decided that catalyst seems like the best approach over raw
# dogging some Modeling Tool Kit. Saves some lines for sure.
# I have no idea if this is right but I believe so.
@variables t [unit = u"s"]
D = Differential(t)
@parameters begin

  pₗᵗᵒᵗ=11.0,          [description="Total lac promoter", connect=Flow, unit=u"nM"]
  pₜᵗᵒᵗ=11.0,         [description="Total Tet Promoter", connect=Flow, unit=u"nM"]
  h₀ = 10.4,         [description="basepair per right hand turn",  unit=u"bp*turn^-1"]
  lₗ=40,              [description="Length of plac", unit=u"bp"]
  lₜ=44,              [description="Length of ptet", unit=u"bp"]
  lₛ=240,             [description="Length of mSpinach and T500 terminator", unit=u"bp"]
  lₘ=63,              [description="Length of MG and T500 terminator", unit=u"bp"]
  lₚ=2892,            [description="Length of plasmid", unit=u"bp"]
  lᵢ = 150,           [description="Length of intergenic region", unit=u"bp"]
  kσₘₘ = 50e-3,       [description="MM constant for supercoiling hillfunctions",unit=u"nM"]
  kₒₚₑₙ=0.04,          [description="Rate of open complex formation", unit=u"s^-1"]
  kᵣ=1/170,           [description="RFP maturation rate", unit=u"s^-1"]
  kaₗ=6e3,             [description="Rate of DNA-free apolacI IPTG binding", unit=u"nM^-1*s^-1"]
  kuaₗ=1,              [description="Rate of apolacI IPTG disassociation", unit =u"s^-1"]
  kbindₗ=10,           [description="lacI-promoter asossiation rate", unit=u"s^-1"]
  kuₗ=0.022,           [description="lacI-promoter disassociation rate", unit=u"s^-1"]
  kaₜ=6e3,             [description="aTc-TetR association rate", unit=u"nM^-1*s^-1"]
  kuaₜ=1,              [description="aTc-TetR disassociation rate", unit=u"s^-1"]
  kbindₜ=10,           [description="tetR-DNA association rate", unit=u"s^-1"]
  kuₜ=0.022,           [description="tetR-DNA disassociation rate", unit=u"s^-1"]
  ρₗ=0,                [description="Rate of lacI production", unit=u"nM*s^-1"]
  ρₜ=0,                [description="Rate of tetR production", unit=u"nM*s^-1"]
  δₛ=log(2)/(30*60),   [description="mSpinach degredation rate", unit=u"s^-1"]
  δₘ=log(2)/(60*60),   [description="MG degredation rate", unit=u"s^-1"]
  δₚ=0,                [description="Average protein degredation rate", unit=u"s^-1"]
  σ₀=-0.065,           [description="Natural B-form DNA supercoil state", unit=u"turn*bp^-1"] 
  σspₗ = σ₀*lₚ/lₗ,       [description="Approximate Optimal supercoiling density, plac"]
  σstₛ = σ₀*lₚ/lₛ,      [description="Approximate Optimal supercoiling density, plac"]
  σspₜ = σ₀*lₚ/lₜ,      [description="Approximate Optimal supercoiling density, pTet"]
  σstₘ = σ₀*lₚ/lₘ,     [description="Approximate Optimal supercoiling density, pTet"]
  fudge = 1,           [description="Fudge Factor", unit=u"nM"]
  gyr₀=18.93, [description="Concentration Gyrase", unit=u"μM"]
  topo₀=2, [description="Conc Topoisomerase", unit=u"μM"]
  τ=0.5, [description="Rate of topoisomerase activity", unit=u"turn*s^-1"]
  γ=0.5, [description="Rate of Gyrase activit", unit=u"turn*s^-1"]
  kσₘₘ=200, [description="Michaelis-Menten constant for gyrase", unit=u"μM"]
end
@variables(begin
  t, [unit=u"s"]
  σtₛ(t), [description="Supercoiling State of mSpinach"], 
  σtₘ(t), [description="Supercoiling State of MG"],
  σpₗ(t), [description="Supercoiling density of plac"],
  σpₜ(t), [description="Supercoiling density of pTet"],
  pₗ(t), [description="conc plac", unit=u"nM"],
  pₗc(t), [description="conc plac-lacI complex", unit=u"nM"],
  pₜ(t), [description="conc pTet", unit=u"nM"],
  pₜc(t), [description="conc pTet-TetR complex", unit=u"nM"],
  Cₛ(t), [description="conc mSpinach", unit=u"nM"],
  Cₘ(t), [description="conc MG", unit=u"nM"],
  ECₛ(t), [description="conc Open Complex for mSpinach", unit=u"nM"],
  CCₛ(t), [description="conc Closed Complex for mSpinach", unit=u"nM"],
  ECₘ(t), [description="conc Open Complex for MG", unit=u"nM"],
  CCₘ(t), [description="conc Closed Complex for MG", unit=u"nM"],
  σ₊(t), [description="Strictly positive compoenent of superoil"], 
  σ₋(t), [description="Strictly negative compoenent of superoil"], 
  Δₖᵢₙₖ(t), [description="kink formed from super coils.", unit=u"bp"],
  nfₛ(t), [description="Length between plac and the kink", unit=u"bp"],
  nfₘ(t), [description="Length between pTet and the kink",unit=u"bp"] 
end)

function m(σ)
  σ₊ = (σ+abs(σ))/2
  σ₋ = (σ-abs(σ))/2
  m=topo₀*τ*(σ₋)/kσₘₘ/(σ₀+abs(σ-σ₀))+gyr₀*γ*(σ₊)/kσₘₘ/(σ₀+abs(σ-σ₀))
end
@register m(σ)

function kᵢₙᵢₜ(σp,σsp)
  kᵢₙᵢₜₘₐₓ = 7e-2
  kᵢₙᵢₜ = σsp*kᵢₙᵢₜₘₐₓ/(σsp+((σp-σsp)^2))
end

function kₑ(σt, σst)
  kₑₘₐₓ = 7e-2
  kₑ = σst*kₑₘₐₓ/(σst +((σt-σst)^2))
end

eqns = [Δₖᵢₙₖ ~ (σtₛ+σtₘ)*h₀,
        nfₛ ~ (lₗ+lₛ+lᵢ/2(pₜ/(pₜ+pₜc))+(lᵢ/2+lₘ)*pₜc/(pₜ+pₜc)-Δₖᵢₙₖ), 
        nfₘ ~ (lₜ+lₘ+lᵢ/2(pₗ/(pₗ-pₗc))+(lᵢ/2+lₛ)*pₗc/(pₗ+pₗc)-Δₖᵢₙₖ),
        D(σtₛ) ~ -(D(Cₛ)-δₛ*Cₛ-D(ECₛ))*(lₛ)/(2*h₀*nfₛ)-(D(ECₛ)-D(CCₛ))*(lₗ/2*h₀*nfₛ)+m(σtₛ),
        D(σpₗ) ~ -(D(ECₛ)-D(CCₛ))*(lₗ/2*h₀*nfₛ)+m(σpₗ),
        D(σtₘ) ~ -(D(Cₘ)-δₛ*Cₘ-D(ECₘ))*(lₘ)/(2*h₀*nfₘ)-(D(ECₘ)-D(CCₘ))*(lₜ/2*h₀*nfₘ)+m(σtₘ),
        D(σpₜ) ~ -(D(ECₘ)-D(CCₘ))*(lₘ/2*h₀*nfₘ)+m(σpₜ)]

@named odesys = ODESystem(eqns,t)

rxn = @reaction_network begin
  ρₗ, 0-->LacI
  δₚ, LacI --> 0
  kₐₗ, + LacI +IPTG --> aLacI
  kᵤₐₗ, aLacI --> LacI+IPTG
  kₛₗ, pₗ + LacI --> pₗC
  kᵤₗ, pₗC --> LacI +pₗ
  kᵢₙᵢₜ(σ, σpsₗ), R + pₗ --> CCₛ
  kₒₚₑₙ, CCₛ --> ECₛ
  kₑ(σ, σtsₛ), ECₛ --> Cₛ+R+pₗ
  δₛ, cₛ --> 0
  ρₜ, 0 --> TetR
  δₚ, TetR --> 0
  kₐₜ, TetR + aTc --> aTetR 
  kᵤₐₜ, aTetR --> aTc + TetR
  kₛₜ, pₜ + TetR --> pₜC 
  kₐₜ, pₜC + aTc --> pₜ+aTetR
  kᵢₙᵢₜ(σ,σpmₜ), R + pₜ --> CCₘ
  kᵣ, CCₘ --> pₜ + R 
  kₒₚₑₙ, CCₘ --> ECₘ
  kₑ(σ, σtmₘ), ECₘ --> Cₘ + R +pₜ
  δₘ, Cₘ --> 0
end
convert(ODESystem, rxn)