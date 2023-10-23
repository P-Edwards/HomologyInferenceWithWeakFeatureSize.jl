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