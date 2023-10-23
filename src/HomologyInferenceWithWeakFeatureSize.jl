module HomologyInferenceWithWeakFeatureSize


# export check_if_vector_is_real
# export check_if_solution_is_degenerate
# export convert_solution_point_to_standard_form
# export apply_distance_squared_to_solution
# export filter_solution_to_bottleneck
# export find_points_on_bottleneck_correspondence
export compute_weak_feature_size
# export parse_bertini_file
# export butterfly_bertini_convert_solution_point_to_standard_form
# export bertini_convert_solution_point_to_standard_form
# export quartic_bertini_convert_solution_point_to_standard_form
export homology_inference
# export minimum_distance_start_system
# export check_if_solution_has_real_endpoints
export subsample_with_function

include("wfs_core.jl")
include("wfs_utilities.jl")
include("sampling_core.jl")
include("utilities_for_examples.jl")

end # module
