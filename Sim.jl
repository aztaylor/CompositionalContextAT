using IJulia, ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
             
@variables begin
  t
  repₗ(t)=0         # conc. of reppressor lacI
  arepₗ(t)=0        # conc. of apolacI
  indᵢ(t)=11       # conc. of inducer IPTG
  promₗ(t)=0        # conc. of lac promoter
  cpromₗ(t)=0       # conc. IPTG bound plac
  ccₛ(t)=0          # lenggth of the closed mSpinach complex.
  ecₛ(t)=0          # length of the mSpinach alongation complex.
  reporterₛ(t)=0    # conc. of the mSplnach reporter
  rnap(t)=18.931   # conc of RNAP
  ribo(t) = 0      # conc of ribosomes
  repₜ(t)=0         # conc of TetR
  indₐ(t)=11       # conc of aTc
  arepₜ(t)=0        # conc of apo TetR
  promₜ(t)=0 
  cpromₜ(t)=0 
  ccₘ(t)=0 
  ecₘ(t)=0 
  σtₛ(t)=-6 
  σtₘ(t)=-3 
  σpₛ(t)=-6 
  σpₘ(t)=-3
end
@parameters begin 
  h₀=10.5 
  lₛ=203 
  lₘ=68 
  lₛ=150 
  lₗ=40 
  lₜ=44 
  lₚ=2892 
  kσₘₘ=50 
  gyr₀=12 
  topo₀=2 
  kₗ₊=7e-2 
  kᵣ=550 
  kₜₓₚ=85 
  kₜₓ=kₜₓₚ/105.5 
  kₒₚₑₙ=0.04 
  kₗₑₐₖ=0.02 
  ρₗ=0 
  ρₜ=0 
  kₜₗ=21/(714+675) 
  kᵪ=1/90 
  kᵣ=1/170 
  k₁₆=0.75 
  kbgl=1.25 
  δₛ=log(2)/(30*60) 
  δₘ=log(2)/(60*60)
  τ=0.5 
  γ=0.5 
  σ₀=-0.065 
  kgyrₘₘ=200 
  kaₗ=6e3 kuaₗ=1 
  kbindₗ=10 
  kuₗ=0.022 
  kaₜ=6e3 
  kuaₜ=1 
  kbindₜ=10 
  kuₜ=0.022 
  δₚ=0
end

xdot=Differential(t)




@begin
