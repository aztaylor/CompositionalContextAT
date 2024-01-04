module UnitfulBio

using Unitful

@unit bp "bp" basePairs 340e-15*u"m" true
@unit nt "nt" nucleotides 340e-15*u"m" true
#@unit turn "turn" Turns 2*π*u"rad" true

Unitful.register(UnitfulBio)
end