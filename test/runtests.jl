using DGGS
using Test

@testset "DGGS.jl" begin
    meta = Dict(
        # specify the operation
        "dggrid_operation" => "GENERATE_GRID",

        # specify the DGG
        "dggs_type" => "FULLER4T",
        "dggs_res_spec" => 0,

        # control the generation
        "clip_subset_type" => "WHOLE_EARTH",
        "geodetic_densify" => 0.0,

        # specify the output
        "cell_output_type" => "KML",
        "cell_output_file_name" => "out",
        "kml_default_color" => "ffff0000",
        "kml_default_width" => 6,
        "kml_name" => "Spherical Icosahedron",
        "kml_description" => "www.discreteglobalgrids.org",
        "kml_default_width" => 6,
        "densification" => 0,
        "precision" => 6
    )

    call_dggrid(meta)
end