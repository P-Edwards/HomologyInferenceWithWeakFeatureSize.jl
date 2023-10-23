using LinearAlgebra, HomotopyContinuation

DEFAULT_THRESHOLD = 1e-8

function total_degree_pattern(F,variables)
    correspondence_polynomials = []
    correspondence_parameters = []
    target_parameters = []
    patch_polynomials = []
    for f in F 
        current_dense_poly,current_dense_params = dense_poly(variables,degree(f))
        current_target_parameters = coeffs_as_dense_poly(f,variables,degree(f))
        random_patch_coefficients = rand(Complex{Float64},length(current_dense_params))
        push!(patch_polynomials,dot(random_patch_coefficients,current_dense_params) - dot(random_patch_coefficients,current_target_parameters))
        push!(target_parameters,current_target_parameters)
        push!(correspondence_polynomials,current_dense_poly)
        push!(correspondence_parameters,current_dense_params)               
    end
    return([ModelKit.Expression(poly) for poly in correspondence_polynomials],[ModelKit.Expression(poly) for poly in patch_polynomials],
        vcat(correspondence_parameters...),vcat(target_parameters...))
end

function check_if_vector_is_real(vec; threshold=DEFAULT_THRESHOLD)
    for coordinate in vec
        if abs(imag(coordinate)) >= threshold
            return false
        end
    end
    return true
end

function check_if_solution_has_real_endpoints(solution; threshold=DEFAULT_THRESHOLD)
    for endpoint in solution["endpoints"]
        if !check_if_vector_is_real(endpoint;threshold=threshold)
            return false
        end
    end
    return true
end

# Convert solution to mulitipliers, endpoints, and circumcenter
function convert_solution_point_to_standard_form(solution,number_of_variables,degree_of_bottlenecks)
    base_for_multipliers = number_of_variables*degree_of_bottlenecks
    coordinates_in_solution = solution[1:base_for_multipliers]
    multipliers = solution[(base_for_multipliers+1):(base_for_multipliers+degree_of_bottlenecks)]
    endpoints = []
    for i in 1:degree_of_bottlenecks
        the_index = i
        arithmetic_progression = (0:(number_of_variables-1))
        arithmetic_progression = [entry*degree_of_bottlenecks + the_index for entry in arithmetic_progression]
        this_endpoint = [coordinates_in_solution[entry] for entry in arithmetic_progression] 
        push!(endpoints,this_endpoint)
    end
    circumcenter = sum([multipliers[i]*endpoints[i] for i in 1:degree_of_bottlenecks])     
    return Dict("endpoints"=>endpoints,"multipliers" => multipliers, "circumcenter" => circumcenter,"distance_squared"=>solution[end])
end

# Check if a solution point is degenerate
function check_if_solution_is_degenerate(solution;threshold=DEFAULT_THRESHOLD)
    for multiplier in solution["multipliers"]
        if abs(multiplier) < threshold
            return true
        end
    end
    endpoints = solution["endpoints"]    
    degree_of_bottlenecks = length(endpoints)
    affine_basis = [endpoints[i]-endpoints[1] for i in 2:degree_of_bottlenecks]
    endpoints_matrix = reduce(hcat,affine_basis)'    
    return rank(endpoints_matrix,atol=threshold) < degree_of_bottlenecks-1
end

function apply_distance_squared_to_solution(solution) 
    if haskey(solution,"distance_squared")
        return solution["distance_squared"]
    end
    first_endpoint = solution["endpoints"][1]
    difference_vector = first_endpoint - solution["circumcenter"]
    return (difference_vector)' * (difference_vector)
end

function real_and_positive_d_value_check(point;threshold=DEFAULT_THRESHOLD)
    return !(abs(imag(apply_distance_squared_to_solution(point))) > threshold) && real(apply_distance_squared_to_solution(point)) > threshold
end