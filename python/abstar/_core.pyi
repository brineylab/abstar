class CoreError(Exception):
    """args == (code: str, message: str)"""

def version() -> str: ...
def fields_toml() -> str: ...
def germline_build(recipe: str, output_dir: str, overwrite: bool, database_paths: dict[str, str]) -> dict: ...
def germline_summary(path: str) -> dict: ...
def germline_gap(queries: list[tuple[str, str]], template_paths: list[str]) -> list[dict]: ...

GAP_DEFAULT_TEMPLATES: list[str]
def run_to_files(
    inputs: list[str],
    records: list[tuple[str, str]],
    output_dir: str,
    db_path: str,
    formats: list[str],
    threads: int | None,
    strict: bool,
    loci: list[str] | None,
    include_nonfunctional: bool,
    sequence_type: str,
) -> dict: ...
def run_to_stdout(
    inputs: list[str],
    records: list[tuple[str, str]],
    db_path: str,
    threads: int | None,
    strict: bool,
    loci: list[str] | None,
    include_nonfunctional: bool,
    sequence_type: str,
) -> dict: ...
def run_in_memory(
    inputs: list[str],
    records: list[tuple[str, str]],
    db_path: str,
    threads: int | None,
    strict: bool,
    loci: list[str] | None,
    include_nonfunctional: bool,
    sequence_type: str,
) -> object: ...
def truth_toml() -> str: ...
def model_summary(path: str) -> dict: ...
def generate_to_files(
    db_path: str, model_path: str, params: dict, output_dir: str, name: str, formats: list[str]
) -> dict: ...
def generate_fasta_to_stdout(db_path: str, model_path: str, params: dict) -> dict: ...
def generate_in_memory(db_path: str, model_path: str, params: dict) -> object: ...
