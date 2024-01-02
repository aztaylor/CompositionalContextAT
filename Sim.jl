using IJulia, ModelingToolkit, DifferentialEquations, Plots, Unitful

@parameters t
D = Differential(t)

@mtkmodel rates begin
  @parameters begin
      lₗ=40, [description="Length of plac", unit=u"bp"]
      lₜ=44, [description="Length of ptet", unit=u"bp"]
      lₛ=240, [description="Length of mSpinach and T500 terminator", unit=u"bp"]
      lₘ=63,  [description="Length of MG and T500 terminator", unit=u"bp"]
      lₚ=2892, [description="Length of plasmid", unit=u"bp"]
      σ₀= -0.065, [description="Natural B-form DNA supercoil state", unit=u"rad*bp^-1"]
      kinitₗ = 7e-2
      kinitₘ = 7e-2, [description="Max initiation rate", unit=u"nt*s^-1"]
      kelongₘ = 7e-2, [description="Max elomgation rate", unit=u"nt*s^-1"]
      kelongₗ = 7e-2
      kσₘₘ= 50, [description="MM constant for supercoiling hillfunctions",unit=u"μM"]
      δₛ=log(2)/(30*60), [description = "mSpinach degredation rate", unit= u"s^-1"]
      δₘ=log(2)/(60*60), [description = "MG degredation rate", unit= u"s^-1"]
      δₚ=0,[description="Average protein degredation rate", unit= u"s^-1"]
      kₒₚₑₙ=0.04, [description = "Rate of open complex formation", unit= u"2*π*rad*s^-1"]
      kᵣ=1/170, [description = "RFP maturation rate", unit= u"s^-1"]
      kaₗ=6e3, [description="Rate of DNA-free apolacI IPTG binding", unit= u"s^-1"]
      kuaₗ=1,[description="Rate of apolacI IPTG disassociation", unit = u"s^-1"]
      kbindₗ=10, [description="lacI-promoter asossiation rate", unit= u"s^-1"]
      kuₗ=0.022, [description ="lacI-promoter disassociation rate", unit= u"s^-1"]
      kaₜ=6e3, [description = "aTc-TetR association rate", unit= u"s^-1"]
      kuaₜ=1, [description="aTc-TetR disassociation rate", unit= u"s^-1"]
      kbindₜ=10, [description="tetR-DNA association rate", unit= u"s^-1"]
      kuₜ=0.022, [description = "tetR-DNA disassociation rate", unit= u"s^-1"]
      ρₗ=0, [description="Rate of lacI production", unit= u"s^-1"]
      ρₜ=0, [description="Rate of tetR production", unit= u"s^-1"] 
      σspₗ = σ₀*lₚ/lₗ, [description="Approximate Optimal supercoiling density, plac"]
      σstₗ = σ₀*lₚ/lₛ, [description="Approximate Optimal supercoiling density, plac"]
      σspₘ = σ₀*lₚ/lₜ, [description="Approximate Optimal supercoiling density, pTet"]
      σstₘ = σ₀*lₚ/lₘ, [description="Approximate Optimal supercoiling density, pTet"]
  end
  @variables begin
    σtₛ(t)=-6, [description="supercoil state of mSpinach ORF"]
    σtₘ(t)=-3, [description="supercoil state of MG ORF"]
    σpₗ(t)=-6, [description="supercoil state of mSpinach promoter"]
    σpₜ(t)=-3, [description="supercoil state of MG promoter"]
  end
  @equations begin
      kinitₗ ~ σspₗ*kinitₘ/(1+(σpₗ-σspₗ)^2)
      kelongₗ ~ σstₗ*kelongₘ/(1+(σtₛ-σstₗ)^2)
      kinitₘ ~ σspₘ*kinitₘ/(1+(σpₜ-σspₘ)^2)
      kelongₘ ~ σstₘ*kelongₘ/(1+(σtₘ-σstₘ)^2)
  end
end

@mtkmodel dnaComplexDynamics begin
  @extend rates()
  @parameters begin
    rnapᵗ=0.0, [description="Total RNAP, in myTXTL no degredataion", connect = Flow]
    promₗᵗ=0.0, [description="Total lac promoter", connect = Flow]
    promₗᵗ=0.0, [description="Total Tet Promoter", connect = Flow]
    reprₗᵗ=0.0, [description="Total Lac Repressor, in myTXTL no degredataion", connect = Flow]
    reprₜᵗ=0.0, [description="Total Tet Repressor, in myTXTL no degredataion", connect = Flow]
    indᵢᵗ=0.0, [description="Total IPTG, not metabolized by the reaction volume", connect = Flow]
    indₐᵗ=0.0, [description="Total aTc, not metabolized by the reaction volume", connect = Flow]
  end
  @variables begin
    ecₛ(t)=0, [description="Number of mSpinach elongation comlexes"]
    ecₘ(t)=0, [description="Number of mSpinach closed dna comlexes"]
    ccₛ(t)=0, [description="Number of MG elongation comlexes"]
    ccₘ(t)=0, [description="Number of MG closed dna comlexes"]
    cpromₗ(t)=0, [description="conc plac-lacI complex" ,unit=u"nM"]
    cpromₜ(t)=0, [description="conc pTet-TetR complex", unit=u"nM"]
    reprₗ(t)=0, [description="conc LacI Repressor", unit=u"nM"]
    reprₜ(t)=0, [description="conc TetR Repressor", unit=u"nM"]   
    areprₗ(t)=0, [description="conc apo LacI", unit=u"nM"]
    areprₜ(t)=0, [description="conc apo TetR", unit=u"nM"]
    indᵢ(t)=indᵢᵗ, [description="conc IPTG", unit=u"nM"]
    indₐ(t)=indₐᵗ, [description="conc aTc", unit=u"ng/μl"]
  end
  @equations begin
    rnapₗᵗ ~ rnap+ecₛ+ecₘ+ccₛ+ccₘ
    promₗᵗ ~ promₗᵗ+ccₛ+ecₛcpromₗ
    promₜᵗ ~ promₜ+ccₘ+ecₘ+cpromₜ
    indᵢᵗ ~ indᵢ+areprₗ+cpromₗ
    indₐᵗ ~ indₐ+areprₜ+cpromₜ
  end
end 

@mtkmodel reporterDynamics begin
  @extend dnaComplexDynamics()
  @variables begin
    reporterₛ(t)=0, [description="mSpinach Transcript", unit="μM"]
    reporterₘ(t)=0, [description="MG Transcript", unit="μM"]
  end
  @equations begin
    D(reporterₛ)~kinitₗ*ecₛ-δₛ*reporterₛ
    D(reporterₘ)~kinitₘ*ecₘ-δₘ*reporterₘ
    D(ecₛ)~kₒₚₑₙ*ccₛ-kinitₗ*ecₛ
    D(ecₘ)~kₒₚₑₙ*ccₘ-kinitₘ*ecₘ
    D(ccₛ)~kelongₗ*(rnapᵗ-ecₛ-ecₘ-ccₛ-ccₘ)*(pₗᵗ-ccₛ-ecₛ-cpromₗ)-(kelongₗ-kₒₚₑₙ)*ccₛ
    D(ccₘ)~kelongₘ*(rnapᵗ-ecₘ-ecₛ-ccₛ-ccₘ)*(pₘᵗ-ccₘ-ecₘ-cpromₜ)-(kelongₘ-kₒₚₑₙ)*ccₘ
    D(reprₗ)~ρₗ+kuaₗ*(indᵢᵗ-indᵢ)+kuₗ(reprₗᵗ-reprₗ-indᵢᵗ-indᵢ)-kaₗ*reprₗ*indᵢ-kbindₗ*reprₗ-σₚ*reprₗ
    D(reprₜ)~ρₜ+kuaₜ*(indₐᵗ-indₐ)*kuₜ(reprₜᵗ-reprₜ-indₐᵗ-indₐ)-kaₜ*reprₜ*ind-kbindₜ*reprₜ-σₚ*reprₜ
    D(indᵢ)~kaₗ*(reprₗ+cpromₗ)*indᵢ+kuaₗ(reprₗᵗ-reprₗ-cpromₗ)
    D(indₐ)~kaₜ*(reprₜ+cpromₜ)*indₐ+kuaₜ(reprₜᵗ-reprₜ-cpromₜ)
  end
end

@mtkmodel nᵢ begin
  @extend reporterDynamics()
  @parameters begin
    lᵢ=203, [description="Lessngth of intergenic region" unit=u"bp"]
    lₛ=240, [description="Length of mSpinach and T500 terminator", unit=u"bp"]
    lₘ=63,  [description="Length of MG and T500 terminator", unit=u"bp"]
    lₗ=40, [description="Length of plac", unit=u"bp"]
    lₜ=44, [description="Length of ptet", unit=u"bp"] 
    Δₖᵢₙₖ = (σtₛ+σtₘ)*h₀, [description="Position of the kink serperating coils from either side", unit=u"bp"] 
  end
  @equations begin
    nfₛ = max((lₗ+lₛ+lᵢ/2(promₜ/(promₜ+cpromₜ))+(lᵢ/2+lₘ)*cpromₜ/(promₜ+cpromₜ)-Δₖᵢₙₖ),0) 
    nfₘ = max((lₜ+lₘ+lᵢ/2(promₗ/(promₗ-cpromₗ))+(lᵢ/2+lₛ)*cpromₗ/(promₗ+cpromₗ)-Δₖᵢₙₖ),0)
  end
end

@mtkmodel σMaintenanceDynamics begin
  @extend reporterDynamics()
  @parameters begin
    gyr₀=12, [description="Concentration Gyrase", unit=u"μM"]
    topo₀=2, [description="Conc Topoisomerase", unit=u"μM"]
    τ=0.5, [description="Rate of topoisomerase activity", unit=u"s^-1"]
    γ=0.5, [description="Rate of Gyrase activit", unit=u"s^-1"]
    kgyrₘₘ=200, [description="Michaelis-Menten constant for gyrase", unit=u"μM"]
    σ₀=-0.065, [description="standard supercoil state"]#, unit=ub"bp"
  end 
  @parameters begin
    σpₗ₊ = 0
    σpₗ₋ = 0
  end  
  @equations begin
    if σpₗ > 0
      σpₗ₊ = σpₗ₊+σpₗ
    elseif σpₗ < 0 
      σpₗ₋ = σpₗ₋+σpₗ
    end
    σtₗ₊ = 0
    σtₗ₋ = 0 
    if σtₛ > 0
       σtₗ₊ = σtₗ₊+σtₛ
    elseif σtₛ < 0
      σtₗ₋ = σtₗ₋+σtₛ; end
    σpₜ₊ = 0
    σpₜ₋ = 0
    if σpₜ > 0
       σpₜ₊ = σpₜ₊+σpₜ
    elseif σpₜ < 0
       σpₜ₋ = σpₜ₋+σpₜ
    end
    σtₘ₊ = 0
    σtₘ₋ = 0
    if σtₘ > 0
       σtₘ₊ = σtₘ₊+σtₘ
    elseif σtₘ < 0
       σtₘ₋ = σtₘ₋+σtₘ
    end
      
    mpₗ~topo₀*τ*(σpₗ₋)/(σ₀+(σpₗ-σ₀)^2)*gyr₀*γ*(σpₗ₊)/(σ₀+(σpₗ-σ₀)^2)
    mtₛ~topo₀*τ*(σtₛ₋)/(σ₀+(σtₛ-σ₀)^2)*gyr₀*γ*(σtₛ₊)/(σ₀+(σtₛ-σ₀)^2)
    mpₜ~topo₀*τ*(σpₜ₋)/(σ₀+(σpₜ-σ₀)^2)*gyr₀*γ*(σpₜ₊)/(σ₀+(σpₜ-σ₀)^2)
    mtₘ~topo₀*τ*(σtₘ₋)/(σ₀+(σtₘ-σ₀)^2)*gyr₀*γ*(σtₘ₊)/(σ₀+(σtₘ-σ₀)^2)
  end
end

@mtkmodel σDynamics begin
  @extend σMaintenanceDynamics()
  @equations begin
    D(σtₛ) ~ -(reporterₛ-δₛ*reporterₛ-D(ecₛ))*(lₛ)/(2*h₀*nfₛ)...
            -(D(ecₛ)-D(ccₛ))*(lₗ/2*h₀*nfₛ)+mpₛ
    D(σpₗ) ~ -(D(ecₛ)-D(ccₛ))*(lₗ/2*h₀*nfₛ)+mpᵢ
    D(σpₜ) ~ (reporterₘ-δₛ*reporterₘ-D(ecₘ))*(lₘ)/(2*h₀*nfₘ)...
            -(D(ecₘ)-D(ccₘ))*(lₗ/2*h₀*nfₘ)+mpₘ
    D(σtₘ) ~ -(D(ecₘ)-D(ccₘ))*(lₗ/2*h₀*nₘ)+mpₜ
  end
end
@parameters begin 
  kₗ₊=7e-2 # Rate forward transcription
  kᵣ=550 # Reverse transcription rate
  kₜₓₚ=85 # Perbase transcription rate
  kₜₓ=kₜₓₚ/105.5 # Average transcription rate
  kₗₑₐₖ=0.02 # Rate of terminator escaping transcription.
end


sys = σDynamics
@mtkbuild model = sys()
#@named sys = ODESystem(eqs,)

