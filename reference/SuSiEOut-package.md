# SuSiEOut: fine-mapping when causal variants may be omitted

SuSiEOut implements SuSiE-OUT, a summary-statistic fine-mapping model
for loci at which a causal variant may have been removed from the tested
set, for example by a minor allele frequency (MAF) threshold or a
SNP-only filter. Each sparse effect can be null, assigned to a tested
variant, or assigned to an omitted variant (an *OUT state*). An omitted
variant is represented by its LD fingerprint: its correlations with the
tested variants, measured in a sequencing reference panel that contains
it.

## Main functions

- [`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md)
  fits the model at one locus.

- [`susieout_fit_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit_grid.md)
  and
  [`susieout_fit_nuisance_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit_nuisance_grid.md)
  average fits over polygenic and LD-mismatch settings.

- [`susieout_candidate_support()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_candidate_support.md)
  and
  [`susieout_credible_sets()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_credible_sets.md)
  summarize a fit.

- [`susieout_simulate_locus()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_simulate_locus.md)
  generates loci with known causal variants.

## Naming

Arguments and fields named `examined_*` refer to tested variants;
`excluded_*` and `outside_*` refer to omitted variants.

## Acknowledgements

SuSiE-OUT builds directly on the Sum of Single Effects (SuSiE) model of
Wang, Sarkar, Carbonetto and Stephens (2020) and its summary-statistic
extension SuSiE-RSS (Zou, Carbonetto, Wang and Stephens, 2022), and it
reuses `susieR` as a reference implementation. We are grateful to the
SuSiE authors for making their methods and software openly available.
SuSiEOut is an independent project and is not affiliated with or
endorsed by the SuSiE authors or the `susieR` maintainers.

## References

Wang G, Sarkar A, Carbonetto P, Stephens M (2020). A simple new approach
to variable selection in regression, with application to genetic fine
mapping. *Journal of the Royal Statistical Society: Series B* 82,
1273–1300. [doi:10.1111/rssb.12388](https://doi.org/10.1111/rssb.12388)

Zou Y, Carbonetto P, Wang G, Stephens M (2022). Fine-mapping from
summary data with the "Sum of Single Effects" model. *PLOS Genetics* 18,
e1010299.
[doi:10.1371/journal.pgen.1010299](https://doi.org/10.1371/journal.pgen.1010299)

## See also

Useful links:

- <https://github.com/weinstocklab/SuSiEOut>

- <https://weinstocklab.github.io/SuSiEOut/>

- Report bugs at <https://github.com/weinstocklab/SuSiEOut/issues>

## Author

**Maintainer**: Josh Weinstock <josh.weinstock@emory.edu>

Authors:

- Josh Weinstock <josh.weinstock@emory.edu>
