using IJulia, ModelingToolkit, DifferentialEquations, Plots, Unitful

@variables begin
  ribo(t) = 0 # conc of ribosomes
end

@parameters begin 
  kₗ₊=7e-2 # Rate forward transcription
  kᵣ=550 # Reverse transcription rate
  kₜₓₚ=85 # Perbase transcription rate
  kₜₓ=kₜₓₚ/105.5 # Average transcription rate
  
  kₗₑₐₖ=0.02 # Rate of terminator escaping transcription.

  kₜₗ=21/(714+675) # RFP/CFP average translation rate
  kᵪ=1/90 # CFP maturation rate

  k₁₆=0.75 # BCD16 binding rate (relative)
  kbgl=1.25 # BCDbgl binding rate (relative)

  σ₀=-0.065 # standard supercoil state
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
      cpromₗ(t)=0, [description="conc plac-lacI complex" ,unit=u"nM"]
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
  @variables begin
    σpₗ(t)
    σtₛ(t)
    σpₜ(t)
    σtₘ(t)
  end

  @parameters begin
      lₗ=40, [description="Length of plac"] #, unit=ub"bp"]
      lₜ=44, [description="Length of ptet"]
      lₚ=2892, [description="Length of plasmid"]#, unit=ub"bp"]
      σ₀= -0.065, [description="Length of plac"]#, unit=ub"bp"]
      kinitₘ = 7e-2, [description="Max initiation rate"]#, unit=u"nt/s"]
      kelongₘ = 7e-2, [description="Max elomgation rate"]#, unit=u"nt/s"]
      kσₘₘ= 50, [description="MM constant for supercoiling hillfunctions",unit=u"μM"]
      σspₗ = σ₀*lₚ/lₗ
      σstₗ = σ₀*lₚ/lₛ
      σspₘ = σ₀*lₚ/lₜ
      σstₘ = σ₀*lₚ/lₘ
  end

  @equations begin
      kinitₗ ~ σspₗ*kinitₘ/(1+(σpₗ(t)-σspₗ)^2)
      kelongₗ ~ σstₗ*kelongₘ/(1+(σtₛ(t)-σstₗ)^2)
      kinitₘ ~ σspₘ*kinitₘ/(1+(σpₜ(t)-σspₘ)^2)
      kelongₘ ~ σstₘ*kelongₘ/(1+(σtₘ(t)-σstₘ)^2)
  end
end

@mtkmodel reporterDynamics begin
  @components begin
      r = rates()
      cLaws = conservationLaws(t)
  end
  @variables begin
      reporterₛ(t)=0, [description="mSpinach Transcript"]
      reporterₘ(t)=0, [description="MG Transcript"]
  end
  @parameters begin
    δₛ=log(2)/(30*60), [description = "mSpinach degredation rate"]
    δₘ=log(2)/(60*60), [description = "MG degredation rate"]
    δₚ=0,[description="Average protein degredation rate"]
    kₒₚₑₙ=0.04, [description = "Rate of open complex formation"]
    kᵣ=1/170, [description = "RFP maturation rate"]
    kaₗ=6e3, [description="Rate of DNA-free apolacI IPTG binding"]
    kuaₗ=1,[description="Rate of apolacI IPTG disassociation"]
    kbindₗ=10, [description="lacI-promoter asossiation rate"]
    kuₗ=0.022, [description ="lacI-promoter disassociation rate"]
    kaₜ=6e3, [description = "aTc-TetR association rate"]
    kuaₜ=1, [description="aTc-TetR disassociation rate"]
    kbindₜ=10, [description="tetR-DNA association rate"]
    kuₜ=0.022, [description = "tetR-DNA disassociation rate"]
    ρₗ=0, [description="Rate of lacI production"]
    ρₜ=0, [description="Rate of tetR production"]    
  end

  @equations begin
    D(reporterₛ)~r.kinitₗ*cLaws.ecₛ-δₛ*reporterₛ
    D(reporterₘ)~r.kinitₘ*cLaws.ecₘ-δₘ*reporterₘ
    D(cLaws.ecₛ)~kₒₚₑₙ*cLaws.ccₛ-r.kinitₗ*cLaws.ecₛ
    D(cLaws.ecₘ)~kₒₚₑₙ*cLaws.ccₘ-r.kinitₘ*cLaws.ecₘ
    D(cLaws.ccₛ)~r.kelongₗ*(rnapᵗ-cLaws.ecₛ-cLaws.ecₘ-cLaws.ccₛ-cLaws.ccₘ)...
                  *(pₗᵗ-cLaws.ccₛ-cLaws.ecₛ-cpromₗ)-(r.kelongₗ-kₒₚₑₙ)*cLaws.ccₛ
    D(cLaws.ccₘ)~r.kelongₘ*(rnapᵗ-cLaws.ecₘ-cLaws.ecₛ-cLaws.ccₛ-cLaws.ccₘ)...
                  *(pₘᵗ-cLaws.ccₘ-cLaws.ecₘ-cpromₜ)-(r.kelongₘ-kₒₚₑₙ)*cLaws.ccₘ
    D(cLaws.repₗ)~ρₗ+kuaₗ*(indᵢᵗ-indᵢ)+kuₗ(repₗᵗ-repₗ-indᵢᵗ-indᵢ)-kaₗ*repₗ*indᵢ...
                  -kbindₗ*repₗ-σₚ*repₗ
    D(cLaws.repₜ)~ρₜ+kuaₜ*(indₐᵗ-indₐ)*kuₜ(repₜᵗ-repₜ-indₐᵗ-indₐ)-kaₜ*repₜ*indₐ...
                  -kbindₜ*repₜ-σₚ*repₜ
    D(cLaws.indᵢ)~kaₗ*(repₗ+cpromₗ)*indᵢ+kuaₗ(repₗᵗ-repₗ-cpromₗ)
    D(cLaws.indₐ)~kaₜ*(repₜ+cpromₜ)*indₐ+kuaₜ(repₜᵗ-repₜ-cpromₜ)
  end
end

@mtkmodel mGyrTopo begin
  @components begin
    s = σdynamics()
  @parameters begin
      gyr₀=12, [description="Concentration Gyrase", unit=u"μM"]
      topo₀=2, [description="Conc Topoisomerase", unit=u"μM"]
      τ=0.5, [description="Rate of topoisomerase activity", unit=u"s^-1"]
      γ=0.5, [description="Rate of Gyrase activit", unit=u"s^-1"]
      kgyrₘₘ=200, [description="Michaelis-Menten constant for gyrase", unit=u"μM"]
      σ₀=-0.065, [description="standard supercoil state"]#, unit=ub"bp"
  end
      
  @variables begin
      σpₗ(t) = σ₀
      σtₛ(t) = σ₀
      σpₜ(t) = σₒ
      σtₘ(t) = σ₀
  end

  @equations begin
    σpₗ₊ = 0;σpₗ₋ = 0 
    if σpₗ > 0; σpₗ₊ = σpₗ₊+σpₗ
    elseif σpₗ < 0; σpₗ₋ = σpₗ₋+σpₗ; end
    σtₗ₊ = 0;σtₗ₋ = 0 
    if σtₛ > 0; σtₗ₊ = σtₗ₊+σtₛ
    elseif σtₛ < 0; σtₗ₋ = σtₗ₋+σtₛ; end
    σpₜ₊ = 0; σpₜ₋ = 0;
    if σpₜ > 0; σpₜ₊ = σpₜ₊+σpₜ
    elseif σpₜ < 0; σpₜ₋ = σpₜ₋+σpₜ; end
    σtₘ₊ = 0; σtₘ₋ = 0;
    if σtₘ > 0; σtₘ₊ = σtₘ₊+σtₘ
    elseif σtₘ < 0; σtₘ₋ = σtₘ₋+σtₘ; end
      
    mpₛ~topo₀*τ*(σpₗ₋)/(σ₀+(σpₗ-σ₀)^2)*gyr₀*γ*(σpₗ₊)/(σ₀+(σpₗ-σ₀)^2)
    mtₛ
    mpₘ
    mtₘ
  end
end

@mtkmodel σDynamics begin
  @components begin
    exdy = reporterDynamics()#, [description="Expression Dynamics"]
    r = rates()
    cLaws=conservationLaws()
    m = mGyrTopo()c
  end
  @variables begin
    σtₛ(t)=-6, [description="supercoil state of mSpinach ORF"]
    σtₘ(t)=-3, [description="supercoil state of MG ORF"]
    σpₗ(t)=-6, [description="supercoil state of mSpinach promoter"]
    σpₜ(t)=-3, [description="supercoil state of MG promoter"]

  end

  @parameters begin
    h₀=10.5, [description="BP per right hand turn of B-DNA"]
    lᵢ=203, [description="Length of intergenic region"]
    nfₛ = lᵢ+r.lₛ
    nfₘ = lᵢ+r.lₘ

  @equations begin
    D(σtₛ) ~ -(cLaws.reporterₛ-cLaws.δₛ*cLaws.reporterₛ-D(cLaws.ecₛ))*(lₛ)/(2*h₀*nfₛ)...
            -(D(cLaws.ecₛ)-D(cLaws.ccₛ))*(lₗ/2*h₀*nfₛ)+m.mpₛ

  end

    
  end
end
      
 # @equations begin
 #     if σpₗ > 0

 #     mpₗ ~ topo₀
#  end