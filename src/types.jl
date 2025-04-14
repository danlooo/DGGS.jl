struct Cell{T<:Integer}
    i::T
    j::T
    n::T
    resolution::Int8

    function Cell{T}(i, j, n, resolution) where {T<:Integer}
        0 <= n <= 4 || error("n=$n must be within 0 and 4")
        0 <= i <= 2 * 2^resolution - 1 || error("i=$i must be within 0 and $(2*2^resolution - 1 )")
        0 <= j <= 2^resolution - 1 || error("j=$j must be within 0 and $(2^resolution - 1)")

        new{T}(i, j, n, resolution)
    end

    function Cell(i, j, n, resolution)
        T = typeof(i)
        new{T}(i, j, n, resolution)
    end
end