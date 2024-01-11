using Catalyst, ModelingToolkit, Plots, Latexify# , Main.UnitfulBio], # ful], # fulAngles
# Have decided that catalyst seems like the best approach over raw
# dogging some Modeling Tool Kit. Saves some lines for sure.
# I have no idea if this is right but I believe so.
@parameters begin
  pₗᵗᵒᵗ=11.0 #,         [description="Total lac promoter", connect=Flow], # =u"nM"]
  pₜᵗᵒᵗ=11.0 #,         [description="Total Tet Promoter", connect=Flow], # =u"nM"]
  h₀ = 10.4 #,          [description="basepair per right hand turn"], #,  unit=u"bp*turn^-1"]
  lₗ=40 #,              [description="Length of plac"], # =u"bp"]
  lₜ=44 #,              [description="Length of ptet"], # =u"bp"]
  lₛ=240 #,             [description="Length of mSpinach and T500 terminator"], # =u"bp"]
  lₘ=63 #,              [description="Length of MG and T500 terminator"], # =u"bp"]
  lₚ=2892 #,            [description="Length of plasmid"], # =u"bp"]
  lᵢ = 150 #,           [description="Length of intergenic region"], # =u"bp"]
  kσₘₘ = 50e-3 #,       [description="MM constant for supercoiling hillfunctions"], # unit=u"nM"]
  kₒₚₑₙ=0.04 #,          [description="Rate of open complex formation"], # =u"s^-1"]
  kᵣ=1/170 #,           [description="RFP maturation rate"], # =u"s^-1"]
  kaₗ=6e3 #,             [description="Rate of DNA-free apolacI IPTG binding"], # =u"nM^-1*s^-1"]
  kuaₗ=1 #,              [description="Rate of apolacI IPTG disassociation"], #  =u"s^-1"]
  kbindₗ=10 #,           [description="lacI-promoter asossiation rate"], # =u"s^-1"]
  kuₗ=0.022 #,           [description="lacI-promoter disassociation rate"], # =u"s^-1"]
  kaₜ=6e3 #,             [description="aTc-TetR association rate"], # =u"nM^-1*s^-1"]
  kuaₜ=1 #,              [description="aTc-TetR disassociation rate"], # =u"s^-1"]
  kbindₜ=10 #,           [description="tetR-DNA association rate"], # =u"s^-1"]
  kuₜ=0.022 #,           [description="tetR-DNA disassociation rate"], # =u"s^-1"]
  ρₗ=0 #,                [description="Rate of lacI production"], # =u"nM*s^-1"]
  ρₜ=0 #,                [description="Rate of tetR production"], # =u"nM*s^-1"]
  δₛ=log(2)/(30*60) #,   [description="mSpinach degredation rate"], # =u"s^-1"]
  δₘ=log(2)/(60*60) #,   [description="MG degredation rate"], # =u"s^-1"]
  δₚ=0 #,                [description="Average protein degredation rate"], # =u"s^-1"]
  σ₀=-0.065 #,           [description="Natural B-form DNA supercoil state"], # =u"turn*bp^-1"] 
  σspₗ = σ₀*lₚ/lₗ #,       [description="Approximate Optimal supercoiling density, plac"], # =2*π*u"rad"]
  σstₛ = σ₀*lₚ/lₛ #,      [description="Approximate Optimal supercoiling density, plac"], # =2*π*u"rad"]
  σspₜ = σ₀*lₚ/lₜ #,      [description="Approximate Optimal supercoiling density, pTet"], # =2*π*u"rad"]
  σstₘ = σ₀*lₚ/lₘ #,     [description="Approximate Optimal supercoiling density, pTet"], # =2*π*u"rad"]
  fudge = 1 #,           [description="Fudge Factor"], # =u"nM"]
  gyr₀=18.93 #, [description="Concentration Gyrase"], # =u"μM"]
  topo₀=2 #, [description="Conc Topoisomerase"], # =u"μM"]
  τ=0.5 #, [description="Rate of topoisomerase activity"], # =u"s^-1"]
  γ=0.5 #, [description="Rate of Gyrase activit"], # =u"s^-1"]
  kσₘₘ=200 #, [description="Michaelis-Menten constant for gyrase"] # =u"μM"]
end
@variables begin
  t 
  σtₛ(t)=σ₀#, [description="Supercoiling State of mSpinach"], # =2*π*u"rad"], 
  σtₘ(t)=σ₀ #, [description="Supercoiling State of MG"], # =2*π*u"rad"],
  σpₗ(t)=σ₀#, [description="Supercoiling density of plac"], # =2*π*u"rad"],
  σpₜ(t)=σ₀ #, [description="Supercoiling density of pTet"], # =2*π*u"rad"],
  pₗ(t)=11 #, [description="conc plac"], # =u"nM"],
  pₗc(t)=0 #, [description="conc plac-lacI complex"], # =u"nM"],
  pₜ(t)=11#, [description="conc pTet"], # =u"nM"],
  pₜc(t)=0 #, [description="conc pTet-TetR complex"], # =u"nM"],
  Cₛ(t)=0 #, [description="conc mSpinach"], # =u"nM"],
  Cₘ(t)=0 #, [description="conc MG"], # =u"nM"],
  ECₛ(t)=0 #, [description="conc Open Complex for mSpinach"], # =u"nM"],
  CCₛ(t)=0 #, [description="conc Closed Complex for mSpinach"], # =u"nM"],
  ECₘ(t)=0 #, [description="conc Open Complex for MG"], # =u"nM"],
  CCₘ(t)=0 #, [description="conc Closed Complex for MG"], # =u"nM"],
  LacI(t)=0 #, [description="conc Lac Repressor", unit=u"nM"],
  aLacI(t)=0 #, [description="Concentration of apoLacI", unit=u"nM"],
  TetR(t)=0 #, [description="Concentration of Tet Repressor", unit=u"nM"],
  aTetR(t)=0 #, [description="Concentration of apoTetR", unit =u"nM"],
  IPTG(t)=1 #, [description="Concentration of IPTG inducer, unit=u"nM"],
  aTc(t) = 1 #, [description="Concentration of aTc Inducer", unit="nM"],
  R(t) = 100 #, [description="Cocentration of RNA Polymerase", unit=u"nM"],
  σ₊(t)=0 #, [description="Strictly positive compoenent of superoil"], # =2*π*u"rad"], 
  σ₋(t)=0 #, [description="Strictly negative compoenent of superoil"], # unit=2*π*u"rad"], 
  Δₖᵢₙₖ(t)=0 #, [description="kink formed from super coils."], # =u"bp"],
  nfₛ(t)=75 #, [description="Length between plac and the kink"], # =u"bp"],
  nfₘ(t)=75 #, [description="Length between pTet and the kink"] #,unit=u"bp"] 
end

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

D = Differential(t)
eqns = [Δₖᵢₙₖ ~ (σtₛ+σtₘ)*h₀,
        nfₛ ~ (lₗ+lₛ+lᵢ/2(pₜ/(pₜ+pₜc))+((lᵢ/2)+lₘ)*pₜc/(pₜ+pₜc)-Δₖᵢₙₖ), 
        nfₘ ~ (lₜ+lₘ+lᵢ/2(pₗ/(pₗ-pₗc))+((lᵢ/2)+lₛ)*pₗc/(pₗ+pₗc)-Δₖᵢₙₖ),
        D(σtₛ) ~ -(D(Cₛ)-δₛ*Cₛ-D(ECₛ))*(lₛ)/(2*h₀*nfₛ)-(D(ECₛ)-D(CCₛ))*(lₗ/2*h₀*nfₛ)+m(σtₛ),
        D(σpₗ) ~ -(D(ECₛ)-D(CCₛ))*(lₗ/2*h₀*nfₛ)+m(σpₗ),
        D(σtₘ) ~ -(D(Cₘ)-δₛ*Cₘ-D(ECₘ))*(lₘ)/(2*h₀*nfₘ)-(D(ECₘ)-D(CCₘ))*(lₜ/2*h₀*nfₘ)+m(σtₘ),
        D(σpₜ) ~ -(D(ECₘ)-D(CCₘ))*(lₘ/2*h₀*nfₘ)+m(σpₜ)]

@named odesys = ODESystem(eqns, t)

rxn = @reaction_network begin
  #@variables σpₗ(t) σpₜ(t) σₛ(t) σₘ(t)
  (ρₗ, δₚ), ∅ <-->LacI
  (kₐₗ,kᵤₐₗ), LacI +IPTG <--> aLacI
  (kₛₗ, kᵤₗ), pₗ + LacI <--> pₗc + IPTG
  (kᵢₙᵢₜ($σpₗ, $σspₗ), kᵣ), R + pₗ --> CCₛ
  kₒₚₑₙ, CCₛ --> ECₛ
  kₑ($σtₛ, $σstₛ), ECₛ --> Cₛ+R+pₗ
  δₛ, Cₛ --> ∅
  (ρₜ, δₚ), ∅ <--> TetR
  (kₐₜ, kᵤₐₜ), TetR + aTc <--> aTetR 
  (kₛₜ, kᵤₜ), pₜ + TetR <--> pₜc + aTc 
  (kᵢₙᵢₜ($σpₜ, $σspₜ), kᵣ), R + pₜ <--> CCₘ
  kₒₚₑₙ, CCₘ --> ECₘ
  kₑ($σtₘ, σstₘ), ECₘ --> Cₘ + R +pₜ
  δₘ, Cₘ --> ∅
end
@named sys = extend(odesys, rxn)

odeprob = ODEProblem(sys, [], [0,3*60*60])

