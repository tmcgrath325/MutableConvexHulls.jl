const DOWN = (0, -1)

function jarvis_search(current, δ, points)
    sin_angles = [normalized_cross2d(δ, p .- current) for p in points]
    min_idx = argmin(sin_angles)
    @assert sin_angles[min_idx] >= 0
    return points[min_idx]
end

function jarvis_march(points)
    # initialize hull
    hull = typeof(points)()
    # choose first point on the hull
    push!(hull, minimum(points))
    # start with the -y direction
    δ = DOWN
    # perform jarvis march 
    current = first(hull)
    while length(hull) < 1 || current !== first(hull)
        push!(hull, jarvis_search(current, δ, points))
        current = last(hull)
        δ = current .- hull[end-1]
    end
    return hull
end