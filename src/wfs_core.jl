using HomotopyContinuation, DynamicPolynomials, LinearAlgebra, IterTools

function bottleneck_correspondence_equations(F,x,k)
    # F - system 
    # x - array of variables for F
    # k - the degree of the bottlenecks
    @assert k > 1
    N = length(x) 
    n = length(F)
    
    # New variables for individual "fiber product" 
    # copies of original system F
    @polyvar y[1:k,1:N]
    correspondence_original_functions = [subs(F,x => y[i,:]) for i in 1:k]
    correspondence_original_functions = vcat(correspondence_original_functions...)
    # This sets up the barycentric coordinates
    # and putative bottleneck z
    @polyvar t[1:k]
    barycenter_equation = [sum(t) - 1]
    z = sum([t[i]*y[i,:] for i in 1:k])
    
    
    # Projective variables and patches
    @polyvar lambdas[0:n,1:k]
    patch_equations = [dot(rand(n+1)+im*rand(n+1),lambdas[:,i]) - 1 for i in 1:k]
    

    # Normal equations 
    gradient_of_original = differentiate(F,x)
    
    gradients_of_factors = [subs(gradient_of_original,x=>y[i,:]) for i in 1:k]
    unflattened_normals = [lambdas[1,i]*(y[i,:] - z) + sum(lambdas[j+1,i]*gradients_of_factors[i][j,:] for j in 1:n) for i in 1:k]
    flattened_normals = vcat((unflattened_normals...)...)
    

    # Distance equations 
    distance_equations = [dot(y[i,:] - z,y[i,:] - z) - dot(y[i+1,:] - z,y[i+1,:] - z) for i in 1:(k-1)]

    @polyvar D
    squared_distance_to_first_endpoint = [D - dot(y[1,:] - z,y[1,:] - z)]
    
    all_equations = [correspondence_original_functions;barycenter_equation;patch_equations;flattened_normals;distance_equations;squared_distance_to_first_endpoint]    
    return(System(all_equations))
end

function minimum_distance_start_system(F,x)
    N = length(x)
    m = length(F)
    @polyvar lambdas[0:m]
    @polyvar epsilon[1:m]
    @polyvar p[1:N]
    @polyvar D
    perturbed_equations = [F[i] - epsilon[i] for i in 1:m]
    
    gradient_of_original = differentiate(F,x)
    unflattened_normal_equations = [[lambdas[1]*(x[j] - p[j]) for j in 1:N] + sum([[lambdas[i+1]*gradient_of_original[i,j] for j in 1:N] for i in 1:m]) ]    
    patch_equation = [dot(rand(m+1)+im*rand(m+1),lambdas) - 1]
    normal_equations = vcat(unflattened_normal_equations...)
    distance_equation = [D - dot(p-x,p-x)]        
    return System([perturbed_equations;patch_equation;normal_equations;distance_equation],variables=[x;lambdas;D],parameters=[p;epsilon])
end

function check_if_vector_is_real(vec; threshold=1e-7)
    for coordinate in vec
        if abs(imag(coordinate)) >= threshold
            return false
        end
    end
    return true
end

function check_if_solution_has_real_endpoints(solution; threshold=1e-7)
    for endpoint in solution["endpoints"]
        if !check_if_vector_is_real(endpoint;threshold=threshold)
            return false
        end
    end
    return true
end


# Convert solution to mulitipliers, endpoints, and circumcenter
function convert_solution_point_to_standard_form(solution, number_of_variables,degree_of_bottlenecks)
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
function check_if_solution_is_degenerate(solution;threshold=1e-7)
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

function filter_solution_to_bottleneck(solution,number_of_variables,degree_of_bottlenecks,distance_system,start,initial_solution;singular=false,return_points=false,threshold=1e-7,reach=false)
    
    # First series of checks is for real 
    # bottlenecks with positive, real barycentric coordinates
    if check_if_solution_is_degenerate(solution) && !reach
        return(false)
    end

    for multiplier in solution["multipliers"]
        if abs(imag(multiplier)) >= threshold && !singular
            return(false)
        end
    end    

    
    # Check for distance 0 or complex distances, which 
    # we can throw out
    distance_squared = apply_distance_squared_to_solution(solution)
    if abs(imag(distance_squared)) >= threshold
        return(false) 
    end
    if real(distance_squared) < threshold
        return(false)
    end
    combination = solution["circumcenter"]
    # In this case we're on a non-dependent positive dimension component
    # but don't have a real bottleneck, so return distance immediately
    if !check_if_vector_is_real(combination;threshold=threshold)
        if return_points || !singular
            return(false)
        end
        return real(sqrt(distance_squared))
    end
    # Otherwise we can filter out
    # the isolated real algebraic bottlenecks which
    # aren't geometric bottlenecks
    
    # First just toss out entries that aren't actually
    # in the convex hull
    for multiplier in solution["multipliers"]
        if real(multiplier) <= threshold || real(multiplier) > (1-threshold)
            return(false)
        end 
    end
    combination = real(combination)    
    zeroes_for_perturbation = [0.0 for i in 1:(length(parameters(distance_system))-number_of_variables)]
    distance_results = solve(distance_system,initial_solution;start_parameters=start,target_parameters=[combination;zeroes_for_perturbation],tracker_options=TrackerOptions(;parameters=:conservative))
    #distance_solutions = solutions(nonsingular(distance_results))
    distance_solutions = solutions(distance_results)
    distance_solutions = [[dist[1:number_of_variables];dist[end]] for dist in distance_solutions if abs(imag(dist[end])) <= threshold]    
    distance_solutions = [real(dist[end]) for dist in distance_solutions if real(dist[end]) >= threshold]
    distance_to_compare = real(distance_squared)   
    if distance_to_compare - minimum(distance_solutions) > threshold       
        return(false)
    end
    if return_points
        return(combination)
    end
    return(sqrt(minimum(distance_solutions)))
end

function cassini_oval(number_of_points,b,the_vars)
    function cispi(point) 
        return cis(point*pi)
    end
    pts = map(cispi, range(0,2,length=number_of_points+1)[1:end-1])
    pts = [(real(point),imag(point)) for point in pts]
    equations = []
    for point in pts
        push!(equations,sum([(the_vars[j] - point[j])^2 for j in 1:length(the_vars)]))
    end
    final_equation = 1
    for eqn in equations
        final_equation *= eqn
    end
    final_equation = final_equation - b^(2*number_of_points)
    return final_equation
end


function cassini_hypersurface(number_of_points,b,the_vars,rotation=false)
    dimension = length(the_vars)
    pts = [ [0.0 for _ in 1:dimension] for _ in 1:(number_of_points-1) ]
    for i in 1:length(pts)
        pts[i][i] = 1
    end    
    # The last point is the unit vector -(e_1+...+e_n)/sqrt(n)
    the_last_point = [1/sqrt(number_of_points-1) * coordinate for coordinate in -sum(pts)]
    push!(pts,the_last_point)
    # Do a random orthogonal transform to all 
    # the points
    if(rotation!=false) 
        random_matrix = rand(dimension,dimension)
        # Q is now a random orthogonal matrix
        Q,R = qr(random_matrix)
        pts = [Q*point for point in pts]
    end
    equations = []
    for point in pts
        push!(equations,sum([(the_vars[j] - point[j])^2 for j in 1:dimension]))
    end
    final_equation = 1
    for eqn in equations
        final_equation *= eqn
    end
    final_equation = final_equation - b^(2*number_of_points)
    return final_equation
end

function find_points_on_bottleneck_correspondence(F,x,k)
    system_for_correspondence = bottleneck_correspondence_equations(F,x,k)
    solved_bottlenecks = solve(system_for_correspondence, start_system = :polyhedral,tracker_options=TrackerOptions(;parameters=:conservative))
    nonsingular_solutions = solutions(nonsingular(solved_bottlenecks))
    singular_solutions = solutions(singular(solved_bottlenecks))

    number_of_variables = length(x)
    degree_of_bottlenecks = k
    singular_solutions = [convert_solution_point_to_standard_form(point,number_of_variables,degree_of_bottlenecks) for point in singular_solutions]
    nonsingular_solutions = [convert_solution_point_to_standard_form(point,number_of_variables,degree_of_bottlenecks) for point in nonsingular_solutions]
    return Dict("singular"=>singular_solutions,"nonsingular"=>nonsingular_solutions)
end

function real_and_positive_d_value_check(point;threshold=1e-7)
    return !(abs(imag(apply_distance_squared_to_solution(point))) > threshold) && real(apply_distance_squared_to_solution(point)) > threshold
end

function compute_weak_feature_size(F;maximum_bottleneck_order=nothing,threshold=1e-8)
    x = variables(F) 
    if maximum_bottleneck_order == nothing
        maximum_bottleneck_order = length(x) + 1
    end
    
    weak_feature_size = Inf

    distance_system_for_filtering = minimum_distance_start_system(F,x)
    start = randn(ComplexF64,length(x)+length(F))
    #initial_system = subs(distance_system_for_filtering,parameters_for_filtering => start)
    initial_solution = solutions(solve(distance_system_for_filtering,target_parameters=start,start_system = :polyhedral,tracker_options=TrackerOptions(;parameters=:conservative)))

    for k in 2:maximum_bottleneck_order
        computed_points = find_points_on_bottleneck_correspondence(F,x,k)
        singular_points = computed_points["singular"]
        nonsingular_points = computed_points["nonsingular"]
        filtered_nonsingular_points = [point for point in nonsingular_points if check_if_solution_has_real_endpoints(point;threshold=threshold)]

        all_points = [singular_points;filtered_nonsingular_points]
        nondegenerate_points = [point for point in all_points if !check_if_solution_is_degenerate(point;threshold=threshold)]
        real_positive_D_value_points = [point for point in nondegenerate_points if real_and_positive_d_value_check(point;threshold=threshold)]
        can_be_distance_filtered = [point for point in real_positive_D_value_points if check_if_solution_has_real_endpoints(point;threshold=threshold)]
        cannot_be_distance_filtered = [point for point in real_positive_D_value_points if !check_if_solution_has_real_endpoints(point;threshold=threshold)]

        distances_isolated = [Inf]
        if length(can_be_distance_filtered) > 0
            distances_isolated = [filter_solution_to_bottleneck(solution,length(x),k,distance_system_for_filtering,start,initial_solution;threshold=threshold) for solution in can_be_distance_filtered]
            distances_isolated = [dist for dist in distances_isolated if dist!=false]
            if length(distances_isolated) == 0
                distances_isolated = [Inf]
            end
        end
        distances_not_isolated = [Inf]
        if length(cannot_be_distance_filtered) > 0 
            distances_not_isolated = [real(sqrt(apply_distance_squared_to_solution(point))) for point in cannot_be_distance_filtered]
        end      
        weak_feature_size = min(weak_feature_size,minimum(distances_isolated))
        weak_feature_size = min(weak_feature_size,minimum(distances_not_isolated))
    end
    return weak_feature_size
end

function parse_bertini_file(filename;max_norm=1e10)
    base_file = readlines(filename)
    base_file = base_file[3:end]
    parsed_tuples = []
    current_tuple = []
    inf_flag = false
    for line in base_file    
        if line == ""
            if !inf_flag
                push!(parsed_tuples,Tuple(current_tuple))
            end
            current_tuple = []
            inf_flag = false
        else
            both_real_and_imaginary_parts = split(line," ")
            both_real_and_imaginary_parts = [parse(Float64,part) for part in both_real_and_imaginary_parts]
            this_entry = both_real_and_imaginary_parts[1]+im*both_real_and_imaginary_parts[2]
            push!(current_tuple,this_entry)
            if abs(this_entry) >= max_norm
                inf_flag = true 
            end
        end
    end
    return parsed_tuples
end

# Convert solution to mulitipliers, endpoints, and circumcenter
function butterfly_bertini_convert_solution_point_to_standard_form(solution, number_of_variables,degree_of_bottlenecks)
    circumcenter = collect(solution[1:number_of_variables])
    endpoints = [] 
    for i in 2:(degree_of_bottlenecks+1)
        push!(endpoints,collect(solution[((i-1)*number_of_variables+1):i*number_of_variables]))
    end
    multipliers = collect(solution[(end-degree_of_bottlenecks+1)+1:end])
    last_multiplier = 1 - sum(multipliers)
    push!(multipliers,last_multiplier)
    return Dict("endpoints"=>endpoints,"multipliers" => multipliers, "circumcenter" => circumcenter)
end

# Convert solution to mulitipliers, endpoints, and circumcenter
function bertini_convert_solution_point_to_standard_form(solution, number_of_variables,degree_of_bottlenecks)    
    endpoints = [] 
    for i in 1:(degree_of_bottlenecks)
        push!(endpoints,collect(solution[((i-1)*number_of_variables+1):i*number_of_variables]))
    end
    circumcenter = collect(solution[(degree_of_bottlenecks*number_of_variables+1):(degree_of_bottlenecks*number_of_variables+number_of_variables)])
    multipliers = collect(solution[(end-degree_of_bottlenecks+1)+1:end])
    last_multiplier = 1 - sum(multipliers)
    push!(multipliers,last_multiplier)
    return Dict("endpoints"=>endpoints,"multipliers" => multipliers, "circumcenter" => circumcenter)
end

# Convert solution to mulitipliers, endpoints, and circumcenter
function quartic_bertini_convert_solution_point_to_standard_form(solution, number_of_variables,degree_of_bottlenecks)    
    endpoints = [] 
    for i in 1:(degree_of_bottlenecks)
        push!(endpoints,collect(solution[((i-1)*number_of_variables+1):i*number_of_variables]))
    end
    #circumcenter = collect(solution[(degree_of_bottlenecks*number_of_variables+1):(degree_of_bottlenecks*number_of_variables+number_of_variables)])
    multipliers = collect(solution[degree_of_bottlenecks*number_of_variables+1:degree_of_bottlenecks*number_of_variables+degree_of_bottlenecks-1])
    last_multiplier = 1 - sum(multipliers)
    push!(multipliers,last_multiplier)
    circumcenter = sum([endpoints[i]*multipliers[i] for i in 1:length(endpoints)])
    return Dict("endpoints"=>endpoints,"multipliers" => multipliers, "circumcenter" => circumcenter)
end

function reach_bertini_convert_solution_point_to_standard_form(solution, number_of_variables,degree_of_bottlenecks=2)
    circumcenter = collect(solution[1:2])
    endpoints = []
    push!(endpoints,collect(solution[3:4]))
    push!(endpoints,collect(solution[5:6]))
    multipliers = [1/2 for i in 1:degree_of_bottlenecks]        
    return Dict("endpoints"=>endpoints,"multipliers" => multipliers, "circumcenter" => circumcenter)
end

