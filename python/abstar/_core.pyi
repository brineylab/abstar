class CoreError(Exception):
    """args == (code: str, message: str)"""

def version() -> str: ...
def fields_toml() -> str: ...
def germline_build(recipe: str, output_dir: str, overwrite: bool, database_paths: dict[str, str]) -> dict: ...
def germline_summary(path: str) -> dict: ...
def run_to_files(
    inputs: list[str],
    output_dir: str,
    db_path: str,
    formats: list[str],
    threads: int | None,
    strict: bool,
    loci: list[str] | None,
) -> dict: ...
def run_to_stdout(
    inputs: list[str], db_path: str, threads: int | None, strict: bool, loci: list[str] | None
) -> dict: ...
def run_in_memory(
    inputs: list[str],
    records: list[tuple[str, str]],
    db_path: str,
    threads: int | None,
    strict: bool,
    loci: list[str] | None,
) -> object: ...
