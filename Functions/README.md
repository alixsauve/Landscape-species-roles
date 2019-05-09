# Functions

* `mw_adjacencies.R`:
mw_adjacencies() returns an array S x S x c (S = number of species, c = number of habitats).
This function takes for input an interaction list formatted as a dataframe with 5 columns
	- `habitat`, the habitat where the interaction occurs;
	- `lower_taxon`, the name of the taxon from the lower trophic level;
	- `upper_taxon`, the name of the taxon from the upper trophic level;
	- `freq_norm`, the normalised frequency of the interaction;
	- `freq`, the raw frequency of the interaction (estimated with the number of individuals caught for instance).
and a species list formatted as a dataframe with at least a column with species name (`taxon`).
mw_adjacencies() also requires three logical variables: `w` to specify whether the metaweb is quantitative or qualitative, `d`, to specify whether the metaweb is directed or not, and `n` to specify whether normalised interaction frequencies should be used.

* `network_role.R` contains multiple functions:
	- `membership()`returns a matrix b of belonging coefficients, dim(b) = S x H, bij being the belonging coefficient of species i to community j.
	- `c_value()` returns the among-habitat connectivity of each species. This function calls membership().
	- `z_value()`returns the within-habitat weight of a species.
	- `z_value_landscape()` returns returns the mean within-habitat weight of all species at the landscape scale. This functions calls the following function: `z_value()`.
	- `classification()` returns the landscape role of a species given its z-c coordinates.

* `random_metaweb.R` contains multiple functions:
	- `full_random_mw()` returns a random metaweb for all types of interactions using `random_mw()`.
	- `random_mw()` returns a random biparite metaweb based on the assumption that insects are more likely to forage where their most abundant resources are.
	- `my_sample()` is identical to the base function `sample()`, except it returns the input set when it contains a single element.

* `feature_scaling.R` returns x after being scaled between a and b.

* `dietBrayCurtis.R` returns the Bray-Curtis dissimilarity between two quantitative diets (vectors).
* `jaccardDistDiet.R` returns the Jaccard distance between two specific diets (vectors).
But see function `vegdist` from the R package `vegan` for pairwise comparisons at the scale of a community.
