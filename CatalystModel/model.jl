using Catalyst, ModelingToolkit, Plots, Latexify
# Have decided that catalyst seems like the best approach over raw
# dogging some Modeling Tool Kit. Saves some lines for sure.
# I have no idea if this is right but I believe so.
@variables t
D = Differential(t)
@parameters begin
  lₗ=40,    [description="Length of plac", unit=u"bp"]
  lₜ=44,    [description="Length of ptet", unit=u"bp"]
  lₛ=240,   [description="Length of mSpinach and T500 terminator", unit=u"bp"]
  lₘ=63,    [description="Length of MG and T500 terminator", unit=u"bp"]
  lₚ=2892,  [description="Length of plasmid", unit=u"bp"]
  kσₘₘ = 50e-3, [description="MM constant for supercoiling hillfunctions",unit=u"nM"]
  kₒₚₑₙ=0.04,[description = "Rate of open complex formation", unit=u"s^-1"]
  kᵣ=1/170, [description = "RFP maturation rate", unit=u"s^-1"]
  kaₗ=6e3,   [description="Rate of DNA-free apolacI IPTG binding", unit=u"nM^-1*s^-1"]
  kuaₗ=1,    [description="Rate of apolacI IPTG disassociation", unit =u"s^-1"]
  kbindₗ=10, [description="lacI-promoter asossiation rate", unit=u"s^-1"]
  kuₗ=0.022, [description ="lacI-promoter disassociation rate", unit=u"s^-1"]
  kaₜ=6e3,   [description = "aTc-TetR association rate", unit=u"nM^-1*s^-1"]
  kuaₜ=1,    [description="aTc-TetR disassociation rate", unit=u"s^-1"]
  kbindₜ=10, [description="tetR-DNA association rate", unit=u"s^-1"]
  kuₜ=0.022, [description = "tetR-DNA disassociation rate", unit=u"s^-1"]
  ρₗ=0,      [description="Rate of lacI production", unit=u"nM*s^-1"]
  ρₜ=0,      [description="Rate of tetR production", unit=u"nM*s^-1"]
  δₛ=log(2)/(30*60),  [description = "mSpinach degredation rate", unit=u"s^-1"]
  δₘ=log(2)/(60*60),  [description = "MG degredation rate", unit=u"s^-1"]
  δₚ=0,       [description="Average protein degredation rate", unit=u"s^-1"]
  σ₀=-0.065,  [description="Natural B-form DNA supercoil state", unit=u"turn*bp^-1"] 
  σspₗ = σ₀*lₚ/lₗ,  [description="Approximate Optimal supercoiling density, plac"]
  σstₛ = σ₀*lₚ/lₛ,  [description="Approximate Optimal supercoiling density, plac"]
  σspₜ = σ₀*lₚ/lₜ,  [description="Approximate Optimal supercoiling density, pTet"]
  σstₘ = σ₀*lₚ/lₘ, [description="Approximate Optimal supercoiling density, pTet"]
end
@equations begin
  D(σtₛ) ~ -(D(reporterₛ)-δₛ*reporterₛ-D(ecₛ))*(lₛ)/(2*h₀*nfₛ)...
  -(D(ecₛ)-D(ccₛ))*(lₗ/2*h₀*nfₛ)+fudge*mtₛ
  D(σpₗ) ~ -(D(ecₛ)-D(ccₛ))*(lₗ/2*h₀*nfₛ)+fudge*mpₗ
  D(σtₘ) ~ (reporterₘ-δₛ*reporterₘ-D(ecₘ))*(lₘ)/(2*h₀*nfₘ)...
  -(D(ecₘ)-D(ccₘ))*(lₜ/2*h₀*nfₘ)+fudge*mpₘ
  D(σpₜ) ~ -(D(ecₘ)-D(ccₘ))*(lₘ/2*h₀*nₘ)+fudge*mpₜ
end

function kᵢₙᵢₜ(σp,σsp)
  kᵢₙᵢₜₘₐₓ = 7e-2
  kᵢₙᵢₜ = σsp*kᵢₙᵢₜₘₐₓ/(σsp+((σp-σsp)^2))
end

function kₑ(σt, σst)
  kₑₘₐₓ = 7e-2
  kₑ = σst*kₑₘₐₓ/(σst +((σt-σst)^2))
end

rxn = @reaction_network begin
  @variables begin
    σ(t)
    σ(t)
    kᵢₙᵢₜ(σ(t), σpt)
    kₑ(σ(t), σst)
  end
  ρₗ, 0-->LacI
  δₚ, LacI --> 0
  kₐₗ, + LacI +IPTG --> aLacI
  kᵤₐₗ, aLacI --> LacI+IPTG
  kₛₗ, pₗ + LacI --> pₗC
  kᵤₗ, pₗC --> LacI +pₗ
  kᵢₙᵢₜ(σₚₛ, σpsₗ), R + pₗ --> CCₛ
  kₒₚₑₙ, CCₛ --> ECₛ
  kₑ(σₜₛ, σts), ECₛ --> Cₛ+R+pₗ
  δₛ, cₛ --> 0
  ρₜ, 0 --> TetR
  δₚ, TetR --> 0
  kₐₜ, TetR + aTc --> aTetR 
  kᵤₐₜ, aTetR --> aTc + TetR
  kₛₜ, pₜ + TetR --> pₜC 
  kₐₜ, pₜC + aTc --> pₜ+aTetR
  kᵢₙᵢₜ(σₚₘ), R + pₜ --> CCₘ
  kᵣ, CCₘ --> pₜ + R 
  kₒₚₑₙ, CCₘ --> ECₘ
  kₑ(σₜₘ), ECₘ --> Cₘ + R +pₜ
  δₘ, Cₘ --> 0
end