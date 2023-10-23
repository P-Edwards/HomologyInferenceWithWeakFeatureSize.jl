using HomotopyContinuation, DynamicPolynomials, LinearAlgebra, IterTools, Combinatorics


DEFAULT_THRESHOLD = 1e-8

function bottleneck_correspondence_equations(F,x,k;parameters=nothing,extra_equations=[],multi_system=true)
    # F - system or array of k systems
    # x - array of variables for F
    # k - the degree of the bottlenecks
    @assert k > 1
    N = length(x) 
    
    if !multi_system
        systems = [F for i in 1:k]
    else
        @assert length(F) == k
        systems = F
    end
    
    # New variables for individual "fiber product" 
    # copies of original systems in F
    @var y[1:k,1:N]
    correspondence_original_functions = [subs(systems[i],x => y[i,:]) for i in 1:k]
    correspondence_original_functions = vcat(correspondence_original_functions...)
    
    # This sets up the barycentric coordinates
    # and putative bottleneck z
    @var t[1:k]
    barycenter_equation = [sum(t) - 1]
    z = sum([t[i]*y[i,:] for i in 1:k])
    
    
    # Projective variables and patches
    system_index = 1
    patch_equations = []
    max_codimension = maximum(length(system) for system in systems)
    @var lambdas[0:max_codimension,1:k]
    final_lambda_variables = []
    for system in systems
        n = length(system)
        push!(final_lambda_variables,lambdas[1:(n+1),system_index])
        push!(patch_equations,dot(rand(n+1)+im*rand(n+1),lambdas[1:(n+1),system_index]) - 1)
        system_index += 1
    end
    patch_equations = vcat(patch_equations...)
    final_lambda_variables = vcat(final_lambda_variables...)

    # Normal equations 
    system_index = 1
    flattened_normals = []
    for system in systems
        n = length(system)
        gradient_of_original = differentiate(system,x)
        gradients_of_factors = subs(gradient_of_original,x=>y[system_index,:])
        unflattened_normals = lambdas[1,system_index]*(y[system_index,:] - z) + sum(lambdas[j+1,system_index]*gradients_of_factors[j,:] for j in 1:n) 
        push!(flattened_normals,vcat(unflattened_normals...))
        system_index += 1
    end
    flattened_normals = vcat(flattened_normals...)
    
    # Distance equations 
    distance_equations = [dot(y[i,:] - z,y[i,:] - z) - dot(y[i+1,:] - z,y[i+1,:] - z) for i in 1:(k-1)]

    @var D
    squared_distance_to_first_endpoint = [D - dot(y[1,:] - z,y[1,:] - z)]
    
    all_equations = [correspondence_original_functions;barycenter_equation;patch_equations;flattened_normals;distance_equations;squared_distance_to_first_endpoint;extra_equations]    
    return(System(all_equations,variables=[y[:];t;final_lambda_variables;D],parameters=parameters))
end

function minimum_distance_start_system(F,x)
    N = length(x)
    m = length(F)
    @var lambdas[0:m]
    @var epsilon[1:m]
    @var p[1:N]
    @var D
    perturbed_equations = [F[i] - epsilon[i] for i in 1:m]
    
    gradient_of_original = differentiate(F,x)
    unflattened_normal_equations = [[lambdas[1]*(x[j] - p[j]) for j in 1:N] + sum([[lambdas[i+1]*gradient_of_original[i,j] for j in 1:N] for i in 1:m]) ]    
    patch_equation = [dot(rand(m+1)+im*rand(m+1),lambdas) - 1]
    normal_equations = vcat(unflattened_normal_equations...)
    distance_equation = [D - dot(p-x,p-x)]        
    return System([perturbed_equations;patch_equation;normal_equations;distance_equation],variables=[x;lambdas;D],parameters=[p;epsilon])
end


function compute_minimum_distance(distance_system,initial_solution,start,target,number_of_variables;threshold=DEFAULT_THRESHOLD)
    distance_results = solve(distance_system,initial_solution;start_parameters=start,target_parameters=target,tracker_options=TrackerOptions(;parameters=:conservative))
    distance_solutions = solutions(distance_results;only_nonsingular=false)
    distance_solutions = [[dist[1:number_of_variables];dist[end]] for dist in distance_solutions if abs(imag(dist[end])) <= threshold]    
    distance_solutions = [real(dist[end]) for dist in distance_solutions if real(dist[end]) >= threshold]
    return(minimum(distance_solutions))
end

function filter_solution_to_bottleneck(solution,number_of_variables,degree_of_bottlenecks,distance_systems;singular=false,return_points=false,threshold=DEFAULT_THRESHOLD,reach=false)
    
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
    distances_to_compare = []    
    for system in keys(distance_systems)
        distance_system = distance_systems[system]["system"]
        initial_solution = distance_systems[system]["initial"]
        start = distance_systems[system]["start"]
        zeroes_for_perturbation = [0.0 for i in 1:(length(parameters(distance_system))-number_of_variables)]
        target = [combination;zeroes_for_perturbation]  
        push!(distances_to_compare,compute_minimum_distance(distance_system,initial_solution,start,target,number_of_variables;threshold=threshold))
    end
    real_rounded_distance_squared = real(distance_squared)
    for distance_to_compare in distances_to_compare
        if real_rounded_distance_squared - distance_to_compare > threshold       
            return(false)
        end
    end
    if return_points
        return(combination)
    end
    return(sqrt(distances_to_compare[1]))
end

function find_points_on_bottleneck_correspondence(F,x,k;solve_with_monodromy=false,multi_system=true,filtered=true,threshold=DEFAULT_THRESHOLD)
    if !multi_system
        systems = [F for i in 1:k]
    else 
        systems = F
    end
    if !solve_with_monodromy
        system_for_correspondence = bottleneck_correspondence_equations(systems,x,k)
        solved_bottlenecks = solve(system_for_correspondence, start_system = :polyhedral,tracker_options=TrackerOptions(;parameters=:conservative))
    else
        # A remark here: The underlying function ModelKit.dense_poly
        # has a somewhat undocumented and positive feature: 
        # every separate call to it generates unique parameter value
        # names for every parameter. So the following doesn't cause
        # "name collisions" 
        polys,patches,parameters,target = total_degree_pattern(vcat(systems...),x)
        start_index = 0
        split_polys = []
        for system in systems
            start_index += 1
            push!(split_polys,polys[start_index:start_index+length(system)-1])
        end
        system_for_correspondence = bottleneck_correspondence_equations(split_polys,x,k;parameters=parameters,extra_equations=patches)   
        bottleneck_initial_solutions = monodromy_solve(system_for_correspondence) 
        solved_bottlenecks = solve(system_for_correspondence,solutions(bottleneck_initial_solutions);start_parameters=bottleneck_initial_solutions.parameters,target_parameters=target,tracker_options=TrackerOptions(;parameters=:conservative))
    end
    
    nonsingular_solutions = solutions(nonsingular(solved_bottlenecks))
    singular_solutions = solutions(singular(solved_bottlenecks);only_nonsingular=false)

    number_of_variables = length(x)
    degree_of_bottlenecks = k
    singular_solutions = [convert_solution_point_to_standard_form(point,number_of_variables,degree_of_bottlenecks) for point in singular_solutions]
    nonsingular_solutions = [convert_solution_point_to_standard_form(point,number_of_variables,degree_of_bottlenecks) for point in nonsingular_solutions]
    if !filtered
        return Dict("singular"=>singular_solutions,"nonsingular"=>nonsingular_solutions)
    end
    filtered_nonsingular_points = [point for point in nonsingular_solutions if check_if_solution_has_real_endpoints(point;threshold=threshold)]
    all_points = [singular_solutions;filtered_nonsingular_points]
    nondegenerate_points = [point for point in all_points if !check_if_solution_is_degenerate(point;threshold=threshold)]
    real_positive_D_value_points = [point for point in nondegenerate_points if real_and_positive_d_value_check(point;threshold=threshold)]
    can_be_distance_filtered = [point for point in real_positive_D_value_points if check_if_solution_has_real_endpoints(point;threshold=threshold)]
    cannot_be_distance_filtered = [point for point in real_positive_D_value_points if !check_if_solution_has_real_endpoints(point;threshold=threshold)]
    return(Dict("singular"=>cannot_be_distance_filtered,"nonsingular"=>can_be_distance_filtered))
end

function compute_weak_feature_size(F;maximum_bottleneck_order=nothing,threshold=DEFAULT_THRESHOLD,solve_with_monodromy=false,system_variables=nothing,multi_system=false)    
    if system_variables == nothing
        if !multi_system
            system_variables = variables(F)
        else
            system_variables = variables(F[1])
        end
    end

    number_of_variables = length(system_variables)

    
    if maximum_bottleneck_order == nothing
        maximum_bottleneck_order = number_of_variables + 1
    end
    
    if !multi_system
        systems = [F for i in 1:maximum_bottleneck_order]
    else
        systems = F
    end
    @assert maximum_bottleneck_order == length(systems)

    weak_feature_size = Inf

    # set up distance filtering data
    distance_systems = Dict()
    for system in unique(systems)
        distance_system_for_filtering = minimum_distance_start_system(system,system_variables)
        start = randn(ComplexF64,number_of_variables+length(system))    
        initial_solution = solutions(solve(distance_system_for_filtering,target_parameters=start,start_system = :polyhedral,tracker_options=TrackerOptions(;parameters=:conservative)))
        distance_systems[system] = Dict("system"=>distance_system_for_filtering,"initial"=>initial_solution,"start"=>start)
    end

    for k in 2:maximum_bottleneck_order
        # Unordered combinations, so doesn't duplicate
        function_combinations = unique(collect(combinations(systems,k)))        
        for subsystem in function_combinations
            computed_points = find_points_on_bottleneck_correspondence(subsystem,system_variables,k;solve_with_monodromy=solve_with_monodromy,threshold=threshold)
            
            can_be_distance_filtered = computed_points["nonsingular"]
            cannot_be_distance_filtered = computed_points["singular"]

            distances_isolated = [Inf]
            if length(can_be_distance_filtered) > 0
                subsetted_distance_systems = filter( p -> p[1] in subsystem, distance_systems)
                distances_isolated = [filter_solution_to_bottleneck(solution,number_of_variables,k,subsetted_distance_systems;threshold=threshold) for solution in can_be_distance_filtered]
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
    end
    return weak_feature_size
end



