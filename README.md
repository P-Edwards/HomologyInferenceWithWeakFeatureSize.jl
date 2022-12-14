HomologyInferenceWithWeakFeatureSize.jl
=======================================

Version
-------
1.0.0


Copyright (C) 2022 [Parker
Edwards](https://sites.nd.edu/parker-edwards/)


Installation
-------------
From Julia in package mode (press `]` to enter package mode, backspace to exit): `add https://github.com/P-Edwards/HomologyInferenceWithWeakFeatureSize.jl.git`


Usage
------
To use multiple cores, run `export JULIA_NUM_THREADS=<number_to_use>` before starting julia. Examples are available at [https://github.com/P-Edwards/wfs-and-reach-examples](https://github.com/P-Edwards/wfs-and-reach-examples).


This package primarily exposes two functions: 

	compute_weak_feature_size(F;maximum_bottleneck_order=nothing,threshold=1e-8)

Where 
* `F` is a list of polynomials
* For `maximum_bottleneck_order` >= 2 then the function will stop at the indicated order. Default is the number of variables for `F` + 1. This may become very expensive for order > 3. 
* `threshold` is a numerical threshold to determine the value at which quantities are considered indistinguishable from 0.

The output is a number, which is the best lower bound on the weak feature size computed using geometric bottlenecks up to the provided order. 


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


Depedencies
-----------
This package uses both [sampling_varieties](https://github.com/olivergafvert/sampling_varieties) by Oliver Gafvert and [Ripserer.jl](https://mtsch.github.io/Ripserer.jl/dev/) by Matija ??ufar. 

License
-------
This repository is distributed under GPLv3. 