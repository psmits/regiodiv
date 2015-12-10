Gradients
=========

Questions/theoretical framework
-------------------------------

Latitudinal diversity gradients, "normal" or "reversed."

Higher extinction and origination at tropics.

Hypotheses based on diversification

-  Museum vs Cradle
  -  Preserve versus create diversity.
-  Out of the Tropics / In to the Tropics
  - Diversify in one and then move into other.

These hypotheses have distinct predictions about the relative birth/death
probabilities, and migration/extirpation, between the environmental
types.

Greater recruitment and/or loss in tropical or temperate zones?

Correlation between recruitment and loss?

Also, density dependence within and between regions. The idea being that if
diversity has intrisic limits, are they due to within region limits, between
region limits, or global limits. also, which provinces are "most" limited.


Data
----

Paleozoic brachiopod occurrence information. Sourced from Foote/Miller work.

Split based on geographic provinces. Currently N Tropical, N Temperate,
S Tropical and S Temperate (defined at 22degree). Would prefer better structure
to the provinces so that tropical/temperate can be higher-level binary
predictor. 


Approach
--------

Hierarchical hidden markov model with density dependent (regional) survival and
origination/migration.

Density-dependent effect is considered constant over time while intercept is
hierarchical across time.

All regions affect the DD of all other regions. Each province has all provinces
effect it, independently.

Convert probabilities into rates a la Liow et al 2015 Ecology Letters

rate = -log(1 - prob) / stage length
