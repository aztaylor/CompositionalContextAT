using Catalyst, ModelingToolkit, Plots, Latexify

function kᵢₙᵢₜ(σ) begin
  asp = 
  kᵢₙᵢₜ = σs*kᵢₙᵢₜₘₐₓ/(σsp+((σp-σsp)^2))
end

function kₑ(σ) begin
  kₑ = σst*kₑₘₐₓ/(σst +((σt-σst)^2))
end

@parameters begin
  
end

rxn = @reaction_network begin
  ρₗ, 0-->LacI
  δₚ, LacI --> 0
  kₐₗ, + LacI +IPTG --> aLacI
  kᵤₐₗ, aLacI --> LacI+IPTG
  kₛₗ, pₗ + LacI --> pₗC
  kᵤₗ, pₗC --> LacI +pₗ
  kᵢₙᵢₜ(σₚₛ), R + pₗ --> CCₛ
  kₒₚₑₙ, CCₛ --> ECₛ
  kₑ(σₜₛ), ECₛ --> Cₛ+R+pₗ
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