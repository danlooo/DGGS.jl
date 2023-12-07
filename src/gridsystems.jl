function aggregate_cell_cube(xout, xin; agg_func=filter_null(mean))
    fac = ceil(Int, size(xin, 1) / size(xout, 1))
    for j in axes(xout, 2)
        for i in axes(xout, 1)
            iview = ((i-1)*fac+1):min(size(xin, 1), (i * fac))
            jview = ((j-1)*fac+1):min(size(xin, 2), (j * fac))
            xout[i, j] = agg_func(view(xin, iview, jview))
        end
    end
end

function GridSystem(cell_cube::CellCube)
    pyramid = Vector{CellCube}(undef, cell_cube.level)
    pyramid[cell_cube.level] = cell_cube

    for coarser_level in cell_cube.level-1:-1:1
        coarser_cell_array = mapCube(
            aggregate_cell_cube,
            pyramid[coarser_level+1].data,
            indims=InDims(:q2di_i, :q2di_j),
            outdims=OutDims(
                Dim{:q2di_i}(range(0; step=1, length=2^(coarser_level - 1))),
                Dim{:q2di_j}(range(0; step=1, length=2^(coarser_level - 1)))
            )
        )
        coarser_cell_cube = CellCube(coarser_cell_array, coarser_level)
        pyramid[coarser_level] = coarser_cell_cube
    end

    return GridSystem(pyramid)
end

function Base.show(io::IO, ::MIME"text/plain", dggs::GridSystem)
    println(io, "DGGS GridSystem with $(length(dggs.data)) levels")
    Base.show(io, "text/plain", dggs.data[1].data.axes)
end