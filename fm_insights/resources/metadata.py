import hail as hl
from .generic import bucket


def get_trait_summary(pop="ALL", min_n_pop=1):
    ht = hl.import_table(f"gs://{bucket}/metadata/trait_summary.txt", impute=True)
    if min_n_pop > 1:
        ht2 = ht.group_by("trait").aggregate(n_pop=hl.agg.count())
        ht2 = ht2.filter(ht2.n_pop >= min_n_pop)
        ht = ht.filter(hl.is_defined(ht2[ht.trait]))
    if pop != "ALL":
        ht = ht.filter(ht.cohort == pop)
    return ht


def get_sample_size_dict(pop, neff=False):
    ht = get_trait_summary(pop)
    if neff:
        ht = ht.annotate(
            n_samples=hl.if_else(
                hl.is_missing(ht.n_cases),
                ht.n_samples,
                ht.n_samples * (ht.n_cases / ht.n_samples) * (1 - ht.n_cases / ht.n_samples),
            )
        )
    d = ht.select("trait", "n_samples").to_pandas().set_index("trait").T.to_dict("records")[0]
    return hl.literal(d)


def get_trait_mapping_dict(pop):
    ht = get_trait_summary(pop)
    d = ht.select("trait_cohort", "trait").to_pandas().set_index("trait_cohort").T.to_dict("records")[0]
    return hl.literal(d)


def get_n_pop_dict():
    ht = get_trait_summary()
    d = ht.group_by("trait").aggregate(n_pop=hl.agg.count()).to_pandas().set_index("trait").T.to_dict("records")[0]
    return hl.literal(d)
