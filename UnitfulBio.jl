__precompile__(true)
module UnitfulBio

import Unitful
using Unitful: @unit 
export @ub_str

@unit bp_ub "bp" basePairs 1e-12*Unitful.m false 

include("bioMacro.jl")

const localUnits = copy(Unitful.basefactors)
function __init__()
  merge!(Unitful.basefactors, localUnits)
  Unitful.register(UnitfulBio)
end
end