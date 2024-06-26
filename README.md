HomologyInferenceWithWeakFeatureSize.jl
=======================================

Version
-------
1.2.0


Copyright (C) 2024 [Parker
Edwards](https://parkeredw.com)


Installation
-------------
For global installation from Julia in package mode (press `]` to enter package mode, backspace to exit): `add https://github.com/P-Edwards/HomologyInferenceWithWeakFeatureSize.jl.git`

To install as a standalone local package: Clone this repository to a directory. Then, from Julia in package mode: `activate <path/to/cloned/copy/of/root/directory>` followed by `instantiate`. 


Basic usage
------
To use multiple cores, run `export JULIA_NUM_THREADS=<number_to_use>` before starting Julia. Examples are available at [https://github.com/P-Edwards/wfs-and-reach-examples](https://github.com/P-Edwards/wfs-and-reach-examples).


This package primarily exposes two functions: 

	compute_weak_feature_size(F;maximum_bottleneck_order=nothing,threshold=1e-8,solve_with_monodromy=false,system_variables=variables(F),multi_system=false)

Where 
* `F` is a list of polynomials
* For `maximum_bottleneck_order` >= 2 the function will stop at the indicated order. Default is the number of variables for `F` + 1. This may become very expensive for order > 3. 
* `threshold` is a numerical threshold to determine the value at which quantities are considered indistinguishable from 0.
* `solve_with_monodromy` is a flag that toggles between a polyhedral homotopy solving method and a [monodromy based one](https://www.juliahomotopycontinuation.org/HomotopyContinuation.jl/stable/monodromy/). There are currently no correctness guarantees with the monodromy method, but it is faster and has worked experimentally.
* `system_variables` is a vector of variables to use for the original system. Useful for requiring the use of variables that do not occur in the functions in `F`.
* `multi_system` is a flag which switches between single system mode (the default), and multi-system mode (described below).

The output is a number, which is the best lower bound on the weak feature size computed using geometric bottlenecks up to the provided order. 

**Example**: Compute the weak feature size of an ellipsoid, which occurs at a 2-bottleneck.
	
	using HomologyInferenceWithWeakFeatureSize, HomotopyContinuation
	@var x y z
	F = [x^2 + y^2 + z^2/2 + x*z/7 - 1]
	compute_weak_feature_size(F;maximum_bottleneck_order=2)

Output: `0.9950352569779831`


The second function is
	
	homology_inference(F;wfs=nothing,subsampling_coefficient=1.0,homology_up_to_degree=1,maximum_bottleneck_order=nothing,delta=1e-10,threshold=1e-7,return_sample=false)

This computes a dense sample of the space according to the weak feature size, then finds the Betti numbers of the corresponding algebraic manifold V_R(F), provided V(F) is a smooth complete intersection of codimension d and `F` contains `d` equations. Inputs are: 

* `F` a list of polynomials 
* `wfs` in the default case the weak feature size will be computed using the above function up to the provided `maximum_bottleneck_order`. The user can instead provide a precomputed number. 
* `subsampling_coefficient` controls the trade-off between sample density and subsampling after the initial sampling. The default is reasonable. 
* `homology_up_to_degree` controls the maximum homology degree for homology inference. Homology in degree 2 usually requires substantial memory. Homology in degree 3 is unlikely to complete.
* `delta` and `threshold` are technical inputs. The defaults are reasonable. 
* If `return_sample` is set to `true`, the output will include the dense sample of the algebraic manifold that was computed. 


The output is a list with entry: 
* `output["betti_numbers"]` a list of the computed Betti numbers in ascending order of degree
* `output["wfs"]` the computed weak feature size
* `output["proportion_of_wfs"]` multiplying wfs by this entry gives the initial density of the sample that was computed, before subsampling
* `output["diagrams"]` is a list of persistence diagrams in the format of [PersistenceDiagrams.jl](https://github.com/mtsch/PersistenceDiagrams.jl). The diagrams are in ascending order of degree. 


Multi-system mode
-----------------
The package also supports weak feature size type computations between multiple varieties, which is useful for computing, e.g., the distance between two varieties V(G) and V(H)

	compute_weak_feature_size(array_of_systems;maximum_bottleneck_order=nothing,threshold=1e-8,solve_with_monodromy=false,system_variables=variables(F[1]),multi_system=true)

Where 
* `array_of_systems` is a vector with `maximum_bottleneck_order` elements that are lists of polynomials

Depedencies
-----------
This package uses both [sampling_varieties](https://github.com/olivergafvert/sampling_varieties) by Oliver Gafvert and [Ripserer.jl](https://mtsch.github.io/Ripserer.jl/dev/) by Matija Čufar. 

License
-------
This repository is distributed under GPLv3. 
