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
    kinitₘₐₓ = 7e-2, [description="Max initiation rate", unit=u"s^-1"]
    kelongₘₐₓ = 7e-2, [description="Max elongation rate", unit=u"s^-1"]
    kσₘₘ= 50, [description="MM constant for supercoiling hillfunctions",unit=u"μM"]
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
    δₛ=log(2)/(30*60), [description = "mSpinach degredation rate", unit= u"s^-1"]
    δₘ=log(2)/(60*60), [description = "MG degredation rate", unit= u"s^-1"]
    δₚ=0,[description="Average protein degredation rate", unit= u"s^-1"]
    σ₀= -0.065, [description="Natural B-form DNA supercoil state", unit=u"2π*rad*bp^-1"] 
    σspₗ = σ₀*lₚ/lₗ, [description="Approximate Optimal supercoiling density, plac"]
    σstₛ = σ₀*lₚ/lₛ, [description="Approximate Optimal supercoiling density, plac"]
    σspₜ = σ₀*lₚ/lₜ, [description="Approximate Optimal supercoiling density, pTet"]
    σstₘ = σ₀*lₚ/lₘ, [description="Approximate Optimal supercoiling density, pTet"]

  end
  @variables begin
    σtₛ(t)=-6, [description="supercoil state of mSpinach ORF"]
    σtₘ(t)=-3, [description="supercoil state of MG ORF"]
    σpₗ(t)=-6, [description="supercoil state of mSpinach promoter"]
    σpₜ(t)=-3, [description="supercoil state of MG promoter"]
    kinitₗ(t) = 0.5, [description="transcription initiation rate plac", unit=u"s^-1"]
    kinitₜ(t) = 0.5, [description="transcription initiation rate pTet", unit=u"s^-1"]
    kelongₛ(t) = 0.5, [description="transcription elongation rate mSpinach", unit=u"s^-1"]
    kelongₘ(t) = 0.5, [description="transcription elongation rate MG", unit=u"s^-1"]
  end
  @equations begin
    kinitₗ ~ σspₗ*kinitₘₐₓ/(σspₗ+((σpₗ-σspₗ)^2))
    kelongₛ ~ σstₛ*kinitₘₐₓ/(σstₛ +((σtₛ-σstₛ)^2))
    kinitₜ ~ σspₜ*kinitₘₐₓ/(σspₜ+((σpₜ-σspₜ)^2))
    kelongₘ ~ σstₘ*kinitₘₐₓ/(σstₘ+((σtₘ-σstₘ)^2))
  end
end

@mtkmodel dnaComplexDynamics begin
  @extend rates()
  @parameters begin
    rnapᵗᵒᵗ=0.0, [description="Total RNAP, in myTXTL no degredataion", unit=u"nM", connect = Flow]
    promₗᵗᵒᵗ=0.0, [description="Total lac promoter", connect = Flow, unit=u"nM"]
    promₜᵗᵒᵗ=0.0, [description="Total Tet Promoter", connect = Flow, unit=u"nM"]
    reprₗᵗᵒᵗ=0.0, [description="Total Lac Repressor, in myTXTL no degredataion", connect = Flow, unit=u"nM"]
    reprₜᵗᵒᵗ=0.0, [description="Total Tet Repressor, in myTXTL no degredataion", connect = Flow, unit=u"nM"]
    indᵢᵗᵒᵗ=0.0, [description="Total IPTG, not metabolized by the reaction volume", connect = Flow, unit=u"nM"]
    indₐᵗᵒᵗ=0.0, [description="Total aTc, not metabolized by the reaction volume", connect = Flow, unit=u"nM"]
  end
  @variables begin
    ecₛ(t)=0, [description="Number of mSpinach elongation comlexes", unit=u"nM"]
    ecₘ(t)=0, [description="Number of mSpinach closed dna comlexes", unit=u"nM"]
    ccₛ(t)=0, [description="Number of MG elongation comlexes", unit=u"nM"]
    ccₘ(t)=0, [description="Number of MG closed dna comlexes", unit=u"nM"]
    promₗ(t)=promₗᵗᵒᵗ, [description="conc plac" ,unit=u"nM"]
    promₜ(t)=promₜᵗᵒᵗ,[description="conc pTet" ,unit=u"nM"]
    cpromₗ(t)=0, [description="conc plac-lacI complex" ,unit=u"nM"]
    cpromₜ(t)=0, [description="conc pTet-TetR complex", unit=u"nM"]
    reprₗ(t)=0, [description="conc LacI Repressor", unit=u"nM"]
    reprₜ(t)=0, [description="conc TetR Repressor", unit=u"nM"]   
    areprₗ(t)=0, [description="conc apo LacI", unit=u"nM"]
    areprₜ(t)=0, [description="conc apo TetR", unit=u"nM"]
    indᵢ(t)=indᵢᵗᵒᵗ, [description="conc IPTG", unit=u"nM"]
    indₐ(t)=indₐᵗᵒᵗ, [description="conc aTc", unit=u"nM"]
    rnap(t)=rnapᵗᵒᵗ, [description="conc rnap" ,unit=u"nM"]
  end
  @equations begin
    rnapᵗᵒᵗ ~ rnap+ecₛ+ecₘ+ccₛ+ccₘ
    promₗᵗᵒᵗ ~ promₗ+ccₛ+ecₛ+cpromₗ
    promₜᵗᵒᵗ ~ promₜ+ccₘ+ecₘ+cpromₜ
    indᵢᵗᵒᵗ ~ indᵢ+areprₗ+cpromₗ
    indₐᵗᵒᵗ ~ indₐ+areprₜ+cpromₜ
  end
end 

@mtkmodel reporterDynamics begin
  @extend dnaComplexDynamics()
  @extend rates()
  @variables begin
    reporterₛ(t)=0, [description="mSpinach Transcript", unit=u"nM"]
    reporterₘ(t)=0, [description="MG Transcript", unit=u"nM"]
  end
  @equations begin
    D(reporterₛ)~kinitₗ*ecₛ-δₛ*reporterₛ
    D(reporterₘ)~kinitₜ*ecₘ-δₘ*reporterₘ
    D(ecₛ)~kₒₚₑₙ*ccₛ-kinitₗ*ecₛ
    D(ecₘ)~kₒₚₑₙ*ccₘ-kinitₜ*ecₘ
    D(ccₛ)~kelongₛ*(rnapᵗᵒᵗ-ecₛ-ecₘ-ccₛ-ccₘ)*(promₗᵗᵒᵗ-ccₛ-ecₛ-cpromₗ)-(kelongₛ-kₒₚₑₙ)*ccₛ
    D(ccₘ)~kelongₘ*(rnapᵗᵒᵗ-ecₘ-ecₛ-ccₛ-ccₘ)*(promₜᵗᵒᵗ-ccₘ-ecₘ-cpromₜ)-(kelongₘ-kₒₚₑₙ)*ccₘ
    D(reprₗ)~ρₗ+kuaₗ*(indᵢᵗᵒᵗ-indᵢ)+kuₗ*(reprₗᵗᵒᵗ-reprₗ-indᵢᵗᵒᵗ-indᵢ)-kaₗ*reprₗ*indᵢ-kbindₗ*reprₗ-δₚ*reprₗ
    D(reprₜ)~ρₜ+kuaₜ*(indₐᵗᵒᵗ-indₐ)*kuₜ*(reprₜᵗᵒᵗ-reprₜ-indₐᵗᵒᵗ-indₐ)-kaₜ*reprₜ*indₐ-kbindₜ*reprₜ-δₚ*reprₜ
    D(indᵢ)~kaₗ*(reprₗ+cpromₗ)*indᵢ+kuaₗ*(reprₗᵗᵒᵗ-reprₗ-cpromₗ)
    D(indₐ)~kaₜ*(reprₜ+cpromₜ)*indₐ+kuaₜ*(reprₜᵗᵒᵗ-reprₜ-cpromₜ)
  end
end

@mtkmodel nᵢ begin
  @extend reporterDynamics()
  @parameters begin
    h₀= 10.5, [description="basepairs per right-hand turn", unit=u"bp*rad^-1*2π^-1"]
  end
  @variables begin  
    σtₛ(t) = σ₀
    σtₘ(t) = σ₀
    nfₛ(t)=100, [description="length of promoter, mSpinach ORF, terminator and available intergenic space", unit=u"bp"]
    nfₘ(t)=100, [description="length of promoter, MG ORF, terminator and available intergenic space", unit=u"bp"]
  end
  @equations begin
    Δₖᵢₙₖ = (σtₛ+σtₘ)*h₀
    nfₛ = max((lₗ+lₛ+lᵢ/2(promₜ/(promₜ+cpromₜ))+(lᵢ/2+lₘ)*cpromₜ/(promₜ+cpromₜ)-Δₖᵢₙₖ),0) 
    nfₘ = max((lₜ+lₘ+lᵢ/2(promₗ/(promₗ-cpromₗ))+(lᵢ/2+lₛ)*cpromₗ/(promₗ+cpromₗ)-Δₖᵢₙₖ),0)
  end
end

@mtkmodel σMaintenanceDynamics begin
  @extend reporterDynamics()
  @parameters begin
    gyr₀=18.93, [description="Concentration Gyrase", unit=u"μM"]
    topo₀=2, [description="Conc Topoisomerase", unit=u"μM"]
    τ=0.5, [description="Rate of topoisomerase activity", unit=u"turn*s^-1"]
    γ=0.5, [description="Rate of Gyrase activit", unit=u"turn*s^-1"]
    kgyrₘₘ=200, [description="Michaelis-Menten constant for gyrase", unit=u"μM"]
    σ₀=-0.065, [description="standard supercoil state"]
  end 
  @variables begin
    σpₗ₊(t) = 0, [description="Decomposition of the plac supercoiling state into strictly positive parts"] 
    σpₗ₋(t) = 0, [description="Decomposition of the plac supercoiling state into strictly negative parts"]
    σtₛ₊(t) = 0, [description="Decomposition of the mSpinach supercoiling state into strictly positive parts"]
    σtₛ₋(t) = 0, [description="Decomposition of the mSpinach supercoiling state into strictly negative parts"]
    σpₜ₊(t) = 0, [description="Decomposition of the pTet supercoiling state into strictly positive parts"]
    σpₜ₋(t) = 0, [description="Decomposition of the pTet supercoiling state into strictly negative parts"]
    σtₘ₊(t) = 0, [description="Decomposition of the MG supercoiling state into strictly positive parts"]
    σtₘ₋(t) = 0, [description="Decomposition of the MG supercoiling state into strictly negative parts"]
    mpₗ(t) = 0, [description="Maintenance Dynamics from topoisomerase and Gyrase"]
    mtₛ(t) = 0, [description="Maintenance Dynamics from topoisomerase and Gyrase"]
    mpₜ(t) = 0, [description="Maintenance Dynamics from topoisomerase and Gyrase"]
    mtₘ(t) = 0, [description="Maintenance Dynamics from topoisomerase and Gyrase"]
  end  
  @equations begin
    if σpₗ > 0
      σpₗ₊ = σpₗ₊+σpₗ
    elseif σpₗ < 0 
      σpₗ₋ = σpₗ₋+σpₗ
    end
    if σtₛ > 0
       σtₛ₊ = σtₛ₊+σtₛ
    elseif σtₛ < 0
      σtₛ₋ = σtₛ₋+σtₛ; end
    if σpₜ > 0
       σpₜ₊ = σpₜ₊+σpₜ
    elseif σpₜ < 0
       σpₜ₋ = σpₜ₋+σpₜ
    end
    if σtₘ > 0
       σtₘ₊ = σtₘ₊+σtₘ
    elseif σtₘ < 0
       σtₘ₋ = σtₘ₋+σtₘ
    end
      
    mpₗ~topo₀*τ*(σpₗ₋)/kσₘₘ/(σ₀+(σpₗ-σ₀)^2)\kσₘₘ*gyr₀*γ*(σpₗ₊)/kσₘₘ/(σ₀+(σpₗ-σ₀)^2)/kσₘₘ
    mtₛ~topo₀*τ*(σtₛ₋)/kσₘₘ/(σ₀+(σtₛ-σ₀)^2)/kσₘₘ*gyr₀*γ*(σtₛ₊)/kσₘₘ/(σ₀+(σtₛ-σ₀)^2)/kσₘₘ
    mpₜ~topo₀*τ*(σpₜ₋)/kσₘₘ/(σ₀+(σpₜ-σ₀)^2)*gyr₀*γ*(σpₜ₊)/kσₘₘ/(σ₀+(σpₜ-σ₀)^2)/kσₘₘ
    mtₘ~topo₀*τ*(σtₘ₋)/kσₘₘ/(σ₀+(σtₘ-σ₀)^2)*gyr₀*γ*(σtₘ₊)/kσₘₘ/(σ₀+(σtₘ-σ₀)^2)/kσₘₘ
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

sys = σDynamics
@mtkbuild model = sys()
#@named sys = ODESystem(eqs,)

