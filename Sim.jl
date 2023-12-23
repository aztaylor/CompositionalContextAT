using IJulia, ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
             
@variables begin
  t
  repₗ(t)=0 # conc. of reppressor lacI
  arepₗ(t)=0 # conc. of apolacI
  indᵢ(t)=11 # conc. of inducer IPTG
  promₗ(t)=0 # conc. of lac promoter
  cpromₗ(t)=0 # conc. IPTG bound plac
  ccₛ(t)=0 # lenggth of the closed mSpinach complex.
  ecₛ(t)=0 # length of the mSpinach alongation complex.
  reporterₛ(t)=0 # conc. of the mSplnach reporter
  rnap(t)=18.931 # conc of RNAP
  ribo(t) = 0 # conc of ribosomes
  repₜ(t)=0 # conc of TetR
  indₐ(t)=11 # conc of aTc
  arepₜ(t)=0 # conc of apoTetR
  promₜ(t)=0 # conc ptet
  cpromₜ(t)=0 # conc bound ptet
  ccₘ(t)=0 # length of MG closed complex
  ecₘ(t)=0 # length of MG elongation complex
  σtₛ(t)=-6 # supercoil state of mSpinach ORF
  σtₘ(t)=-3 # supercoil state of MG ORF
  σpₛ(t)=-6 # supercoil state of mSpinach promoter
  σpₘ(t)=-3 # supercoil state of MG promoter
end

@parameters begin 
  h₀=10.5 # BP per right hand turn of B-DNA
  lₛ=203 # Length of intergenic region
  lₘ=68 # Length MS
  lₛ=150 # Length of MS
  
  lₜ=44 # Length of ptet

  
  
  kₗ₊=7e-2 # Rate forward transcription
  kᵣ=550 # Reverse transcription rate
  kₜₓₚ=85 # Perbase transcription rate
  kₜₓ=kₜₓₚ/105.5 # Average transcription rate
  kₒₚₑₙ=0.04 # Rate of open complex formation
  kₗₑₐₖ=0.02 # Rate of terminator escaping transcriptio
  ρₗ=0 # Rate of lacI production
  ρₜ=0 # Rate of tetR production
  kₜₗ=21/(714+675) # RFP/CFP average translation rate
  kᵪ=1/90 # CFP maturation rate
  kᵣ=1/170 # RFP maturation rate
  k₁₆=0.75 # BCD16 binding rate (relative)
  kbgl=1.25 # BCDbgl binding rate (relative)
  δₛ=log(2)/(30*60) # mSpinach degredation rate
  δₘ=log(2)/(60*60) # MG degredation rate
  
  σ₀=-0.065 # standard supercoil state
  
  kaₗ=6e3 # Rate of DNA-free apolacI IPTG binding
  kuaₗ=1 # Rate of apolacI IPTG disassociation
  kbindₗ=10 # lacI-promoter asossiation rate
  kuₗ=0.022 # lacI-promoter disassociation rate
  kaₜ=6e3 # aTc-TetR association rate
  kuaₜ=1 # aTc-TetR disassociation rate
  kbindₜ=10 # tetR-DNA association rate
  kuₜ=0.022 # tetR-DNA disassociation rate
  δₚ=0 # Average protein degredation rate
end

D=Differential(t)

@mtkmodel rates begin
    @parameters begin
        lₗ = 40, [description="Length of plac" units=u"bp"]
        lₚ=2892, [description="Length of plasmid" units=u"bp"]
        σ₀= -0.065, [description="Length of plac" units=u"bp"]
        kinitₘ = 7e-2, [description="Max initiation rate" units=u"nt/s"]
        kelongₘ = 7e-2, [description="Max elomgation rate" units=u"nt/s"]
        kσₘₘ= 50, [description="Michaelis-Menten constant for supercoiling hillfunctions", units=u"uM"]

    @structural_paramaters begin
        σspₗ = σ₀*lₚ/lₗ
        σstₗ = σ₀*lₚ/lₛ
        σspₘ = σ₀*lₚ/lₜ
        σstₘ = σ₀*lₚ/lₘ

    @variables
        σpₗ(t)
        σtₗ(t)
        σpₘ(t)
        σtₘ(t)

    @functions
        kinitₗ ~ σspₗ*kinitₘ/(1+(σpₗ(t)-σspₗ)^2)
        kelongₗ ~ σstₗ*kelongₘ/(1+(σtₗ(t)-σstₗ)^2)
        kinitₘ ~ σspₘ*kinitₘ/(1+(σpₘ(t)-σspₘ)^2)
        kelongₘ ~ σstₘ*kelongₘ/(1+(σtₘ(t)-σstₘ)^2)
    end
end

@mtkmodel mGyrTopo begin
    @paramaters begin
        gyr₀=12, [description="Concentration Gyrase" units=u"uM"]
        topo₀=2 [description:"Conc. Topoisomerase" units=u"uM"]
        τ=0.5, [description="Rate of topoisomerase activity" units=u"s^-1"]
        γ=0.5, ["Rate of Gyrase activity
end" units=u"s^-1"]
        kgyrₘₘ=200, [description="Michaelis-Menten constant for gyrase" units=u"uM"]