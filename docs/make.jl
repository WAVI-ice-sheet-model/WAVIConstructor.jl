using Documenter
using WAVIConstructor

makedocs(
    sitename = "WAVIConstructor.jl Documentation",
    modules = [WAVIConstructor],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ]
)
