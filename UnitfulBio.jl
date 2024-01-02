module UnitfulBio

using Unitful

@unit bp "bp" basePairs 340e-15*u"m" false 
@unit nt "nt" nucleotides 340e-15*u"m" false

Unitful.register(UnitfulBio)
end