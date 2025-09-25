using Revise
using LiveServer
using WAVIConstructor

Revise.revise()
include("make.jl")
servedocs(
    foldername="./docs",
    include_dirs=["./src/"]
)
 