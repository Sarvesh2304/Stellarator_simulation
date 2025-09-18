"""
Installation script for StellaratorPhysics package
"""

println("Installing StellaratorPhysics package...")

# Activate the environment
using Pkg
Pkg.activate(".")

# Install dependencies
println("Installing dependencies...")
Pkg.instantiate()

# Test the installation
println("Testing installation...")
try
    include("test/test_basic.jl")
    println("✓ Installation successful!")
    println("You can now use the package with: using StellaratorPhysics")
catch e
    println("⚠ Installation test failed: $(e)")
    println("Please check the error messages above and try again.")
end
