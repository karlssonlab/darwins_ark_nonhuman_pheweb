
# TODO: standardize on `def lookup_{key}(): return overrides[key]`?

from .utils import PheWebError, load_module_from_filepath
from . import parse_utils

import os, boltons.fileutils
from typing import Optional,Any,Dict,Tuple,List


# All state lives in this dictionary.
# PheWeb should run fine when this is empty.
# Modify with `set_override(key,value)`.
overrides: Dict[str,Any] = {}


replacements = {
    'minimum_maf': 'assoc_min_maf',
    'cache': 'cache_dir',
    'limit_num_variants': 'debugging_limit_num_variants',
    'aliases': 'field_aliases',
}
def set_override(key:str, value:Any) -> None:
    if key in replacements:
        key = replacements[key]
        #raise PheWebError("Use `{}` instead of `{}`.".format(replacements[key], key))
    if key=='download_pheno_sumstats' and value=='secret':
        overrides['secret_download_pheno_sumstats'] = True
    if key=='download_top_hits' and value=='hide':
        overrides['show_download_top_hits_button'] = False
    if key=='download_phenotypes' and value=='hide':
        overrides['show_download_phenotypes_button'] = False
    elif key == 'login':
        if not isinstance(value, dict): raise Exception(value)
        overrides['login'] = True
        overrides['login_GOOGLE_LOGIN_CLIENT_ID'] = value['GOOGLE_LOGIN_CLIENT_ID']
        overrides['login_GOOGLE_LOGIN_CLIENT_SECRET'] = value['GOOGLE_LOGIN_CLIENT_SECRET']
        overrides['login_allowlist'] = [email.lower() for email in value.get('whitelist', [])]
    elif key in ['extra_per_variant_fields', 'extra_per_assoc_fields', 'extra_per_pheno_fields', 'null_values']:
        print("Ignoring config for {} because it's not implemented.".format(key))
        #raise PheWebError("parse_utils customization isn't implemented yet.")
    elif key == 'field_aliases':
        if not isinstance(value, dict): raise PheWebError("field_aliases should be dict, not {!r}".format(value))
        for alias,field_name in value.items():
            if not isinstance(alias, str): raise PheWebError("alias {!r} should be str.".format(alias))
            if field_name not in parse_utils.fields:
                raise PheWebError("The field_name {!r} (for alias {!r}) isn't a real field.  Add it or use one of {}.".format(field_name, alias, list(parse_utils.fields)))
        value = {alias.lower(): field_name.lower() for alias,field_name in value.items()}
        overrides[key] = {**parse_utils.default_field_aliases, **value}
    else:
        overrides[key] = value

def load_overrides_from_file(filepath:str) -> None:
    if not os.path.isfile(filepath): raise PheWebError("Cannot load config from non-existent {!r}".format(filepath))
    try:
        module = load_module_from_filepath("config", filepath)
    except Exception:
        raise PheWebError("Failed running {!r}".format(filepath))
    for key in dir(module):
        if not key.startswith('_'):
            set_override(key, getattr(module, key))


# Utilities
def _mkdir_and_check_readable(dirpath:str) -> bool:
    if not os.path.isdir(dirpath):
        try: boltons.fileutils.mkdir_p(dirpath)
        except PermissionError: return False
    return _is_readable(dirpath)
def _is_readable(filepath:str) -> bool: return os.access(filepath, os.R_OK)
def _check_overrides_type(key:str, value_type:type) -> None:
    if key in overrides and not isinstance(overrides[key], value_type):
        raise PheWebError("configuration for {} must be of type {}, but {!r} is of type {}".format(key, value_type, overrides[key], type(overrides[key])))
def _get_config_optional_str(key:str) -> Optional[str]:
    _check_overrides_type(key, str)
    return overrides.get(key, None)
def _get_config_str(key:str, default:Optional[str] = None) -> str:
    _check_overrides_type(key, str)
    if default is None and key not in overrides:
        raise PheWebError("overrides is missing required config for {!r}".format(key))
    return overrides[key] if default is None else overrides.get(key, default)
def _get_config_optional_int(key:str) -> Optional[int]:
    _check_overrides_type(key, int)
    return overrides.get(key, None)
def _get_config_int(key:str, default:int) -> int:
    _check_overrides_type(key, int)
    return overrides.get(key, default)
def _get_config_optional_float(key:str) -> Optional[float]:
    _check_overrides_type(key, float)
    return overrides.get(key, None)
def _get_config_float(key:str, default:float) -> float:
    _check_overrides_type(key, float)
    return overrides.get(key, default)
def _get_config_bool(key:str, default:bool) -> bool:
    _check_overrides_type(key, bool)
    return overrides.get(key, default)


# Core config
def get_data_dir() -> str:
    if 'PHEWEB_DATADIR' in os.environ: data_dir = os.environ['PHEWEB_DATADIR']
    else: data_dir = _get_config_str('data_dir', os.path.curdir)
    data_dir = os.path.abspath(data_dir)
    if not _mkdir_and_check_readable(data_dir):
        raise PheWebError("PheWeb cannot use data_dir {!r}, because it either can't create it or can't read it.".format(data_dir))
    return data_dir
def get_cache_dir() -> Optional[str]:
    key = 'cache_dir'
    if not overrides.get(key):
        return None
    cache_dir = _get_config_str(key, '~/.pheweb/cache')
    cache_dir = os.path.abspath(os.path.expanduser(cache_dir))
    if not _mkdir_and_check_readable(cache_dir):
        print("Warning: caching is disabled because PheWeb either can't create or can't read {!r}.".format(cache_dir))
        overrides[key] = False
        return None
    return cache_dir


## Debugging config
def is_debug_mode() -> bool: return 'PHEWEB_DEBUG' in os.environ or _get_config_bool('debug', False)
def get_debugging_limit_num_variants() -> Optional[int]: return _get_config_optional_int('debugging_limit_num_variants')
def is_allowed_to_download() -> bool: return not _get_config_bool('disallow_downloads', False)

## Loading config
def get_hg_build_number() -> int:
    ret = _get_config_int('hg_build_number', 19)
    if ret not in [19,38]: raise PheWebError("hg_build_number must be either 19 or 38, not {!r}".format(ret))
    return ret
def get_grch_build_number() -> int:
    hg = get_hg_build_number()
    return hg if hg >= 38 else 18 + hg

## Species config
def get_species() -> str:
    '''Which species this data dir serves. `PHEWEB_SPECIES` env var overrides the
    `species` config key; default is the backward-compatible 'dog'.'''
    from . import species as _species
    if 'PHEWEB_SPECIES' in os.environ:
        ret = os.environ['PHEWEB_SPECIES']
    else:
        ret = _get_config_str('species', _species.DEFAULT_SPECIES)
    if ret not in _species.SPECIES:
        raise PheWebError("species must be one of {}, not {!r}".format(_species.known_species(), ret))
    return ret
def get_species_profile() -> Dict[str,Any]:
    '''The active species profile (see pheweb/species.py). `genes_bed`,
    `ref_fasta_pattern` and `recomb_map` may be overridden per data dir by
    config.py keys.'''
    from . import species as _species
    profile = _species.get_profile(get_species())
    override_bed = _get_config_optional_str('genes_bed')
    if override_bed is not None: profile['genes_bed'] = override_bed
    override_fasta = _get_config_optional_str('ref_fasta_pattern')
    if override_fasta is not None: profile['ref_fasta_pattern'] = override_fasta
    override_recomb = _get_config_optional_str('recomb_map')
    if override_recomb is not None: profile['recomb_map'] = override_recomb
    # Branding is a property of the *study*, not the species: several dog datasets
    # can be served from one codebase and should not all announce themselves as the
    # same site. The species profile supplies a sane default; a data dir overrides.
    override_title = _get_config_optional_str('site_title')
    if override_title is not None: profile['site_title'] = override_title
    override_display = _get_config_optional_str('display_name')
    if override_display is not None: profile['display_name'] = override_display
    return profile

def get_recomb_map_filepath() -> Optional[str]:
    '''Absolute path of the local recombination map, or None if this data dir has
    none (no `recomb_map` in the species profile, or the file is missing).

    The existence check is what the region view keys the recombination track on, so
    a profile naming a map that was never built simply hides the track instead of
    serving 404s to LocusZoom.

    A `recomb_map` that is an absolute path (or starts with `~`) is used as-is, so
    several data dirs on the same machine can share one copy of what is really an
    assembly-level annotation; a bare filename resolves inside the data dir.'''
    filename = get_species_profile().get('recomb_map')
    if not filename: return None
    if filename.startswith('~') or os.path.isabs(filename):
        filepath = os.path.abspath(os.path.expanduser(filename))
    else:
        filepath = os.path.join(get_data_dir(), filename)
    # The .tbi is as essential as the map itself -- without it pysam can't do a
    # region query, so treat a missing index as "no map" rather than failing later.
    if not os.path.exists(filepath) or not os.path.exists(filepath + '.tbi'):
        return None
    return filepath

def get_num_procs(cmd:Optional[str] = None) -> int:
    import multiprocessing
    key = 'num_procs'
    try: return int(overrides[key][cmd])
    except Exception: pass
    try: return int(overrides[key]['*'])
    except Exception: pass
    try: return int(overrides[key])
    except Exception: pass
    n_cpus = multiprocessing.cpu_count()
    return 1 if n_cpus==1 else int(n_cpus * 3/4)



## Parsing config
def get_assoc_min_maf() -> float: return _get_config_float('assoc_min_maf', 0)
def get_field_aliases() -> Dict[str,str]:
    return overrides.get('field_aliases', parse_utils.default_field_aliases)


## Manhattan / top-hits / top-loci config
def get_within_pheno_mask_around_peak() -> int: return _get_config_int('within_pheno_mask_around_peak', 500_000)
def get_between_pheno_mask_around_peak() -> int: return _get_config_int('between_pheno_mask_around_peak', 1_000_000)
def get_manhattan_num_unbinned() -> int: return _get_config_int('manhattan_num_unbinned', 500)
def get_manhattan_peak_max_count() -> int: return _get_config_int('manhattan_peak_max_count', 500)
def get_manhattan_peak_pval_threshold() -> float: return _get_config_float('manhattan_peak_pval_threshold', 1e-6)
def get_manhattan_peak_sprawl_dist() -> int: return _get_config_int('manhattan_peak_sprawl_dist', 200_000)
def get_significance_threshold() -> float:
    '''The genome-wide significance p-value for this data dir.

    Default 5e-8, the human GWAS convention, which is not right for every study:
    it assumes a particular number of independent tests. A dog study on a smaller,
    more structured genome may justify a different threshold, so this is per-data-
    dir (`significance_threshold` in config.py).

    Drives the Manhattan significance line and its tooltip, which peaks get gene
    labels, the region view's significance line, and (unless separately overridden)
    `manhattan_peak_variant_counting_pval_threshold`.'''
    ret = _get_config_float('significance_threshold', 5e-8)
    if not 0 < ret < 1:
        raise PheWebError("significance_threshold must be between 0 and 1, not {!r}".format(ret))
    peak_threshold = get_manhattan_peak_pval_threshold()
    if ret >= peak_threshold:
        # manhattan.py asserts this; catch it here with an explanation, because
        # the assert fires deep inside a load step with no message.
        raise PheWebError(
            "significance_threshold ({!r}) must be stricter (smaller) than "
            "manhattan_peak_pval_threshold ({!r}); peaks have to be detected before "
            "their significant variants can be counted. Lower significance_threshold, "
            "or raise manhattan_peak_pval_threshold.".format(ret, peak_threshold))
    return ret
def get_manhattan_peak_variant_counting_pval_threshold() -> float:
    '''Variants stronger than this are counted into `num_significant_in_peak`.
    Defaults to the data dir's significance threshold so one knob moves both;
    override separately only if you need them to differ.'''
    return _get_config_float('manhattan_peak_variant_counting_pval_threshold', get_significance_threshold())
def get_top_hits_pval_cutoff() -> float: return _get_config_float('top_hits_pval_cutoff', 1e-6)


## Pheno correlation config
def should_show_correlations() -> bool: return _get_config_bool('show_correlations', False)
def get_pheno_correlations_pvalue_threshold() -> float: return _get_config_float('pheno_correlations_pvalue_threshold', 0.05)


## Serving config
def get_lzjs_version() -> str: return _get_config_str('lzjs_version', '0.13.0')
def should_allow_variant_json_cors() -> bool: return _get_config_bool('allow_variant_json_cors', True)
def get_urlprefix() -> str: return _get_config_str('urlprefix', '').rstrip('/')
def get_custom_templates_dir() -> Optional[str]:
    key = 'custom_templates'
    custom_templates_dir = _get_config_str(key, 'custom_templates')
    custom_templates_dir = os.path.abspath(os.path.expanduser(custom_templates_dir))
    if not _is_readable(custom_templates_dir):
        return None
    return custom_templates_dir
def is_login_required() -> bool: return bool(overrides.get('login'))
def get_login_google_id_and_secret() -> Tuple[str,str]:
    if not overrides.get('login'): raise PheWebError("Missing login config")
    return (
        _get_config_str('login_GOOGLE_LOGIN_CLIENT_ID'),
        _get_config_str('login_GOOGLE_LOGIN_CLIENT_SECRET')
    )
def get_login_allowlist() -> List[str]:
    key = 'login_allowlist'
    _check_overrides_type(key, list)
    return overrides.get(key, [])  # type: ignore
def get_secret_key() -> str: return _get_config_str('SECRET_KEY', 'nonsecret key')
def should_show_download_top_hits_button() -> bool: return _get_config_bool('show_download_top_hits_button', True)
def should_show_download_phenotypes_button() -> bool: return _get_config_bool('show_download_phenotypes_button', True)
def is_secret_download_pheno_sumstats() -> bool: return _get_config_bool('secret_download_pheno_sumstats', False)
def get_google_analytics_id() -> Optional[str]: return _get_config_optional_str('GOOGLE_ANALYTICS_TRACKING_ID')
def get_sentry_id() -> Optional[str]: return _get_config_optional_str('SENTRY_DSN')
def should_show_manhattan_filter_button() -> bool: return _get_config_bool('show_manhattan_filter_button', False)
def should_show_manhattan_filter_consequence() -> bool: return _get_config_bool('show_manhattan_filter_consequence', False)
