module Gensys
using LinearAlgebra
include("gensysdt.jl")
include("gensysct.jl")
export gensysdt,
gensysct
end
