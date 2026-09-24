
# Tests for the species profile abstraction (pheweb/species.py + conf accessors).
# Dog is the backward-compatible default; cat is selected by the `species` config key.

import pytest


def _clear_chrom_caches():
    from pheweb import utils
    utils.get_chrom_order_list.cache_clear()
    utils.get_chrom_order.cache_clear()
    utils.get_chrom_aliases.cache_clear()


@pytest.fixture(autouse=True)
def _reset_conf():
    # Isolate each test: clear the species/genes_bed overrides and the lru_caches.
    from pheweb import conf
    for key in ('species', 'genes_bed', 'ref_fasta_pattern', 'recomb_map',
                'significance_threshold', 'manhattan_peak_variant_counting_pval_threshold',
                'manhattan_peak_pval_threshold', 'site_title', 'display_name'):
        conf.overrides.pop(key, None)
    _clear_chrom_caches()
    yield
    for key in ('species', 'genes_bed', 'ref_fasta_pattern', 'recomb_map',
                'significance_threshold', 'manhattan_peak_variant_counting_pval_threshold',
                'manhattan_peak_pval_threshold', 'site_title', 'display_name'):
        conf.overrides.pop(key, None)
    _clear_chrom_caches()


def test_default_is_dog():
    from pheweb import conf, utils
    assert conf.get_species() == 'dog'
    assert utils.get_chrom_order_list() == [str(c) for c in range(1, 39)] + ['39']
    prof = conf.get_species_profile()
    assert prof['ucsc_db'] == 'canFam4'
    assert prof['ucsc_chrom_map'] == {'39': 'X'}


def test_cat_profile():
    from pheweb import conf, utils
    conf.set_override('species', 'cat')
    _clear_chrom_caches()
    assert conf.get_species() == 'cat'
    assert utils.get_chrom_order_list() == [str(c) for c in range(1, 20)]
    prof = conf.get_species_profile()
    assert prof['genes_bed'] == 'cat_genes.bed'
    assert prof['ucsc_db'] == 'hub_6476843_GCF_018350175.1'
    assert prof['ucsc_chrom_map']['2'] == 'A2'
    assert prof['ucsc_chrom_map']['19'] == 'X'


def test_genes_bed_config_override():
    from pheweb import conf
    conf.set_override('species', 'cat')
    conf.set_override('genes_bed', 'my_cat_genes.bed')
    assert conf.get_species_profile()['genes_bed'] == 'my_cat_genes.bed'


def test_recomb_map_in_profiles():
    from pheweb import conf
    assert conf.get_species_profile()['recomb_map'] == 'dog_average_canFam4_recomb.tsv.gz'
    conf.set_override('species', 'cat')
    assert conf.get_species_profile()['recomb_map'] is None


def test_recomb_map_config_override():
    from pheweb import conf
    conf.set_override('recomb_map', 'my_map.tsv.gz')
    assert conf.get_species_profile()['recomb_map'] == 'my_map.tsv.gz'


def test_recomb_map_filepath_none_when_absent(tmp_path, monkeypatch):
    '''A profile naming a map that was never built must resolve to None, so the
    region view hides the track instead of pointing LocusZoom at a 404.'''
    from pheweb import conf
    # PHEWEB_DATADIR beats the config override in get_data_dir(), and the
    # container image sets it, so clear it before pointing at tmp_path.
    monkeypatch.delenv('PHEWEB_DATADIR', raising=False)
    conf.set_override('data_dir', str(tmp_path))
    assert conf.get_recomb_map_filepath() is None

    # present but unindexed is also unusable: pysam can't do a region query
    (tmp_path / 'dog_average_canFam4_recomb.tsv.gz').write_bytes(b'')
    assert conf.get_recomb_map_filepath() is None

    (tmp_path / 'dog_average_canFam4_recomb.tsv.gz.tbi').write_bytes(b'')
    assert conf.get_recomb_map_filepath() == str(tmp_path / 'dog_average_canFam4_recomb.tsv.gz')

    conf.overrides.pop('data_dir', None)


def test_recomb_map_absolute_path_is_used_as_is(tmp_path, monkeypatch):
    '''An absolute `recomb_map` lets several data dirs share one copy of what is
    really an assembly-level annotation, rather than duplicating 157MB each.'''
    from pheweb import conf
    shared = tmp_path / 'shared'
    shared.mkdir()
    (shared / 'map.tsv.gz').write_bytes(b'')
    (shared / 'map.tsv.gz.tbi').write_bytes(b'')

    monkeypatch.delenv('PHEWEB_DATADIR', raising=False)
    conf.set_override('data_dir', str(tmp_path / 'elsewhere'))
    conf.set_override('recomb_map', str(shared / 'map.tsv.gz'))
    assert conf.get_recomb_map_filepath() == str(shared / 'map.tsv.gz')

    conf.overrides.pop('data_dir', None)


def test_significance_threshold_default_and_override():
    from pheweb import conf
    assert conf.get_significance_threshold() == 5e-8
    conf.set_override('significance_threshold', 4e-7)
    assert conf.get_significance_threshold() == 4e-7
    conf.overrides.pop('significance_threshold', None)


def test_variant_counting_threshold_follows_significance_threshold():
    '''One knob should move both, so the Manhattan line and num_significant_in_peak
    cannot silently disagree -- while still allowing an explicit split.'''
    from pheweb import conf
    conf.set_override('significance_threshold', 4e-7)
    assert conf.get_manhattan_peak_variant_counting_pval_threshold() == 4e-7
    conf.set_override('manhattan_peak_variant_counting_pval_threshold', 5e-8)
    assert conf.get_manhattan_peak_variant_counting_pval_threshold() == 5e-8
    assert conf.get_significance_threshold() == 4e-7
    for k in ('significance_threshold', 'manhattan_peak_variant_counting_pval_threshold'):
        conf.overrides.pop(k, None)


def test_significance_threshold_rejects_bad_values():
    '''manhattan.py asserts counting < peak-extending; that assert fires deep in a
    load step with no message, so catch the condition at config time instead.'''
    from pheweb import conf
    from pheweb.utils import PheWebError
    # not a probability
    conf.set_override('significance_threshold', 1.5)
    with pytest.raises(PheWebError):
        conf.get_significance_threshold()
    # looser than manhattan_peak_pval_threshold (default 1e-6)
    conf.set_override('significance_threshold', 1e-5)
    with pytest.raises(PheWebError):
        conf.get_significance_threshold()
    # ...but fine if the peak threshold is raised to match
    conf.set_override('manhattan_peak_pval_threshold', 1e-4)
    assert conf.get_significance_threshold() == 1e-5
    for k in ('significance_threshold', 'manhattan_peak_pval_threshold'):
        conf.overrides.pop(k, None)


def test_unknown_species_raises():
    from pheweb import conf
    from pheweb.utils import PheWebError
    conf.set_override('species', 'ferret')
    with pytest.raises(PheWebError):
        conf.get_species()


def test_branding_is_overridable_per_data_dir():
    '''site_title/display_name describe the study, not the species -- one codebase
    serves several dog datasets and they should not all claim the same name.'''
    from pheweb import conf
    prof = conf.get_species_profile()
    assert prof['site_title'] == 'DAP PheWeb'            # dog profile default
    assert prof['display_name'] == 'Dog Aging Project'
    conf.set_override('site_title', "Darwin's Ark Dog Compulsive Disorder PheWeb")
    conf.set_override('display_name', "Darwin's Ark")
    prof = conf.get_species_profile()
    assert prof['site_title'] == "Darwin's Ark Dog Compulsive Disorder PheWeb"
    assert prof['display_name'] == "Darwin's Ark"
    for k in ('site_title', 'display_name'):
        conf.overrides.pop(k, None)
