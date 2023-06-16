using Makie
using GeoMakie

# Can not use Makie plot recipies, because we need to specify the axis for GeoMakie
# see https://discourse.julialang.org/t/accessing-axis-in-makie-plot-recipes/66006

function plot_geo_cube(geo_cube::YAXArray; latitude_name::String="lat", longitude_name::String="lon")
    latitude_axis = getproperty(geo_cube, Symbol(latitude_name))
    longitude_axis = getproperty(geo_cube, Symbol(longitude_name))
    fig = Figure()
    ga1 = GeoAxis(fig[1, 1]; dest="+proj=wintri", coastlines=true)
    sf = surface!(ga1, longitude_axis, latitude_axis, geo_cube.data; colormap=:viridis, shading=false)
    cb1 = Colorbar(fig[1, 2], sf; label="Value", height=Relative(0.5))
    fig
end

function plot_cell_cube(cell_cube::YAXArray)
    geo_cube = get_geo_cube(cell_cube)
    plot_geo_cube(geo_cube)
end

function plot_grid_system(dggs::GridSystem, resolution::Int=3)
    geo_cube = get_geo_cube(dggs, resolution)
    plot_geo_cube(geo_cube)
end