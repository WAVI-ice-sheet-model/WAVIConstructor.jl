using Documenter
using WAVIConstructor

makedocs(
    sitename = "WAVIConstructor.jl Documentation",
    modules = [WAVIConstructor],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    # Only check that exported/public API functions are documented
    # Internal helper functions with docstrings are optional
    checkdocs = :exports,
)
