using HomotopyContinuation, DynamicPolynomials, LinearAlgebra, IterTools, Ripserer, NearestNeighbors

# This function is (C) Oliver Gafvert, obtained and included under GPLv3
# From sampling_varieties: https://github.com/olivergafvert/sampling_varieties/blob/main/src/sampling_varieties.jl
function sampling_fixed_density(F,epsilon)
    # F - system, n - ambient dimension, x - array of variables of F, epsilon - density
    x = variables(F)
    n = length(x)
    @var y[1:n] p[1:n] gamma[1:n]
    d=length(F) # codimension of variety
    k = n-d # dimension of variety
    @var lambda[1:d] # lagrange multipliers

    gradx = differentiate(F, x)

    δ = epsilon/(2*sqrt(n)) # grid size


    # Compute the bounding box by computing the EDD starting from the center of the largest bottleneck

    q = [1 for _ in 1:n]
    system = [F; map(j -> x[j]-q[j]-dot(lambda, gradx[:, j]), 1:n)]
    result = solve(system, start_system = :polyhedral)


    # Extract farthest point from q to X and use as box length

    critical_points = sort!(map(c -> (norm(c[1:n]-q), c[1:n]), real_solutions(nonsingular(result))), by = a -> a[1])
    b = critical_points[end][1]
    indices = [i for i in -b:δ:b];


    # Compute basic sample

    samples = []
    counter = 0

    for s in IterTools.subsets(1:n, k)
        Ft = System([F; map(i -> x[s[i]]-p[i]-q[s[i]], 1:k)]; parameters=p[1:k])
        p₀ = randn(ComplexF64, k)
        result_p₀ = solve(Ft, target_parameters = p₀)
        S_p₀ = solutions(result_p₀)

        slice_samples = solve(
        Ft,
        S_p₀;
        start_parameters =  p₀,
        target_parameters = [map(j -> indices[p1[j]], 1:k) for p1 in Iterators.product(map(j-> 1:length(indices), s)...)],
        transform_result = (r,p) -> real_solutions(r),
        flatten = true
        )
        samples = vcat(samples, slice_samples)
    end


    # Compute extra sample

    extra_samples = []
    extra_counter = 0

    for l in 1:k-1
        for s in IterTools.subsets(1:n, l)
            Ft = [F; map(i -> x[s[i]]-p[i]-q[s[i]], 1:l)]
            grad = differentiate(Ft, x)
            system = System([Ft; map(j -> x[j]-y[j]-dot(gamma[1:n-k+l], grad[:, j]), 1:n)]; parameters=[y; p[1:l]])
            p₀ = randn(ComplexF64, n+l)
            result_p₀ = solve(system, target_parameters = p₀)
            S_p₀ = solutions(result_p₀)

            slice_samples = solve(
            system,
            S_p₀;
            start_parameters =  p₀,
            target_parameters = [[randn(Float64, n); map(j -> indices[p1[j]], 1:l)] for p1 in Iterators.product(map(j-> 1:length(indices), s)...)],
            transform_result = (r,p) -> real_solutions(r),
            flatten = true
            )
            extra_samples = vcat(extra_samples, slice_samples)
        end
    end

    return vcat(samples, extra_samples)
end

function subsample_with_function(sample,func,norm_func=norm;tree=nothing,has_index=false)
    sample = sort(sample,by=func,rev=true)
    if tree == nothing
        if has_index
            tree = KDTree( reduce(vcat,transpose.([point[1:(end-1)] for point in sample]))' )
        else        
            tree = KDTree( reduce(vcat,transpose.(sample))' )
        end
    end
    should_keep = [true for k in 1:length(sample)]
    for k in 1:length(should_keep)
        if !should_keep[k]
            continue
        end
        the_current_point = sample[k]
        the_current_query_point = the_current_point
        if has_index
            the_current_query_point = the_current_query_point[1:(end-1)]
        end
        current_threshold = func(the_current_point)
        indices_to_remove = inrange(tree,the_current_query_point,current_threshold)
        for index in indices_to_remove
            if index > k
                should_keep[index] = false
            end
        end
    end
    output = [sample[k] for k in 1:length(sample) if should_keep[k]]
    return output
end

function constant_func_curry(the_constant=0) 
    return function(x) return the_constant end
end

function naive_one_proportion(F,target_density,feature_size;subsampling_coefficient=1.0)
    x = variables(F) 
    n = length(x)
    original_sample = sampling_fixed_density(F,target_density*feature_size)
    subsampling_allowance = target_density*feature_size*subsampling_coefficient
    the_sample = subsample_with_function(original_sample,constant_func_curry(subsampling_allowance))
    return (the_sample,original_sample)
end 

function an(n) 
    return sqrt((2*n)/(n+1))
end

function feature_size_to_density_constant(n,feature_size,delta)    
    epsilon = (feature_size/2 - delta*an(n))/(2*(an(n)^2 - 1/2))
    return epsilon
end

function feature_size_to_density_proportion(n,feature_size;delta=1e-10,gamma=1)
    delta_w = delta/feature_size
    numerator = 1/2 - delta_w*an(n)
    denominator = (1+gamma)*(2*an(n)^2 - 1)    
    return numerator/denominator
end

function homology_inference(F;wfs=nothing,subsampling_coefficient=1.0,homology_up_to_degree=1,maximum_bottleneck_order=nothing,delta=1e-10,threshold=1e-7,return_sample=false)
    if wfs == nothing
        wfs = compute_weak_feature_size(F;maximum_bottleneck_order=maximum_bottleneck_order,threshold=threshold)
    end
    ambient_dimension = length(variables(F))
    lambda = feature_size_to_density_proportion(ambient_dimension,wfs;delta=delta)
    subsample, = naive_one_proportion(F,lambda,wfs;subsampling_coefficient=subsampling_coefficient)

    # subsample is now a sample dense enough to perform homology inference
    epsilon = 2*lambda*wfs
    starting_radius_for_homology_inference = 2*epsilon
    ending_radius_for_homology_inference = 2*(2*epsilon*an(ambient_dimension)+delta)
    maximum_radius_required = ending_radius_for_homology_inference*1.01
    persistence_diagrams = ripserer(subsample,threshold=maximum_radius_required,dim_max=homology_up_to_degree)    

    betti_numbers = []
    # Count the points contributing to each homology degree
    for diagram in persistence_diagrams
        betti_number_points_in_this_diagram = [point for point in diagram if point[1] <= starting_radius_for_homology_inference && point[2] >= ending_radius_for_homology_inference]
        push!(betti_numbers,length(betti_number_points_in_this_diagram))
    end

    output = Dict("betti_numbers"=>betti_numbers,"wfs"=>wfs,"proportion_of_wfs"=>lambda,"diagrams"=>persistence_diagrams)
    if return_sample
        output["sample"] = subsample
    end
    return output
end