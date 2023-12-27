using IJulia, ModelingToolkit, DifferentialEquations, Plots, Unitful
@variables begin
  t
  ribo(t) = 0 # conc of ribosomes
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
  
  kₗₑₐₖ=0.02 # Rate of terminator escaping transcriptio
  ρₗ=0 # Rate of lacI production
  ρₜ=0 # Rate of tetR production
  kₜₗ=21/(714+675) # RFP/CFP average translation rate
  kᵪ=1/90 # CFP maturation rate
  kᵣ=1/170 # RFP maturation rate
  k₁₆=0.75 # BCD16 binding rate (relative)
  kbgl=1.25 # BCDbgl binding rate (relative)

  
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
@variables t
D=Differential(t)

@connector conservationLaws begin

  @variables begin
      t
      rnapᵗ(t) = 0.0, [connect=Flow, description="Total RNAP"]
      pₗᵗ(t)=0.0, [connect=Flow, description="Total lac promoter"]
      pₗᵗ(t)=0.0, [connect=Flow, description="T [otal Tet Promoter"]
      repₗᵗ(t)=0.0, [connect=Flow, description="Total Lac Repressor"]
      repₜᵗ(t)=0.0, [connect=Flow, description="Total Tet Repressor"]
      indᵢᵗ(t)=0.0, [connect=Flow, description="Total IPTG"]
      indₐᵗ(t)=0.0, [connect=Flow, description="Total aTc"]
      ccₛ(t)=0, [description="Number of mSpinach elongation comlexes"]
      ecₛ(t)=0, [description="Number of mSpinach closed dna comlexes"]
      ccₘ(t)=0, [description="Number of MG elongation comlexes"]
      ecₘ(t)=0, [description="Number of MG closed dna comlexes"]
      rnap(t)=0, [description="conc RNAP", unit=u"nM"]
      promₗ(t)=0, [description="conc plac", unit=u"nM"]
      cprom(t)=0, [description="conc plac-lacI complex" ,unit=u"nM"]
      promₜ(t)=0, [description="conc pTet", unit=u"nM"]
      cpromₜ(t)=0, [description="conc pTet-TetR complex", unit=u"nM"]
      repₗ(t)=0, [description="conc lacI Repressor", unit=u"nM"] 
      arepₗ(t)=0, [description="conc apo LacI repressor", unit=u"nM"] 
      repₜ(t)=0, [description="conc TetR Repressor", unit=u"nM"]
      arepₜ(t)=0, [description="conc apo TetR", unit=u"nM"]
      indᵢ(t)=0, [description="conc IPTG", unit=u"nM"]
      indₐ(t)=0, [description="conc aTc", unit=u"ng/μl"]
  end
  @equations begin
      rnapₗᵗ ~ rnap+ecₛ+ecₘ+ccₛ+ccₘ
      promₗᵗ ~ promₗᵗ+ccₛ+ecₛcpromₗ
      promₜᵗ ~ promₜ+ccₘ+ecₘ+cpromₜ
      indᵢᵗ ~ indᵢ+arepₗ+cpromₗ
      indₐᵗ ~ indₐ+arepₜ+cpromₜ
  end
end

@mtkmodel rates begin
  @parameters begin
      lₗ = 40, [description="Length of plac", unit=u"bp"]
      lₚ=2892, [description="Length of plasmid", unit=u"bp"]
      σ₀= -0.065, [description="Length of plac", unit=u"bp"]
      kinitₘ = 7e-2, [description="Max initiation rate", unit=u"nt/s"]
      kelongₘ = 7e-2, [description="Max elomgation rate", unit=u"nt/s"]
      kσₘₘ= 50, [description="Michaelis-Menten constant for supercoiling hillfunctions", unit=u"uM"]
      σspₗ = σ₀*lₚ/lₗ
      σstₗ = σ₀*lₚ/lₛ
      σspₘ = σ₀*lₚ/lₜ
      σstₘ = σ₀*lₚ/lₘ
  end

  @variables begin
      σpₗ(t)
      σtₗ(t)
      σpₘ(t)
      σtₘ(t)
  end

  @functions begin
      kinitₗ ~ σspₗ*kinitₘ/(1+(σpₗ(t)-σspₗ)^2)
      kelongₗ ~ σstₗ*kelongₘ/(1+(σtₗ(t)-σstₗ)^2)
      kinitₘ ~ σspₘ*kinitₘ/(1+(σpₘ(t)-σspₘ)^2)
      kelongₘ ~ σstₘ*kelongₘ/(1+(σtₘ(t)-σstₘ)^2)
  end
end

@mktmodel reporterDynamics begin
  @components begin
      r = rates()
      cLaws = conservationLaws()
  end
  @variables begin
      reporterₛ(t)=0
      reporterₘ(t)=0
      cLaws.ecₛ(t)=0
      cLaws.ecₘ(t)=0
      cLaws.ccₛ(t)=0
      cLaws.ccₘ(t)=0
      cLaws.repₗ(t)=0
      cLaws.repₜ(t)=0
      cLaws.indᵢ(t)=1
      cLaws.indₐ(t)=100
  end
  @parameters begin
      δₛ=log(2)/(30*60) # mSpinach degredation rate
      δₘ=log(2)/(60*60) # MG degredation rate
      kₒₚₑₙ=0.04 # Rate of open complex formation
      r.kii
      
  end
end

@mtkmodel σDynamics begin
  @variables begin
  σtₛ(t)=-6 # supercoil state of mSpinach ORF
  σtₘ(t)=-3 # supercoil state of MG ORF
  σpₛ(t)=-6 # supercoil state of mSpinach promoter
  σpₘ(t)=-3 # supercoil state of MG promoter   
  end
end

@mtkmodel mGyrTopo begin
  @paramaters begin
      gyr₀=12, [description="Concentration Gyrase", unit=u"uM"]
      topo₀=2, [description:"Conc. Topoisomerase", unit=u"uM"]
      τ=0.5, [description="Rate of topoisomerase activity", unit=u"s^-1"]
      γ=0.5, [description="Rate of Gyrase activit", unit=u"1/s"]
      kgyrₘₘ=200, [description="Michaelis-Menten constant for gyrase", unit=u"uM"]
      σ₀=-0.065, [definition="standard supercoil state", unit=u"bp"]
  end
      
  @variables begin
      σpₗ(t)
      σtₗ(t)
      σpₘ(t)
      σtₘ(t)
      
  end
      
  @functions begin
      if σpₗ > 0

      mpₗ ~ topo₀
  end