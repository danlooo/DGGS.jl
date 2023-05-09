using DGGS
using Test

@testset "DGGS.jl" begin
    using DGGRID7_jll
    dggrid() do dggrid_path
        run(`env`)
        run(`$dggrid_path`)
    end
end

const dggrid_binpath = @generate_wrappers(DGGRID7_jll)
