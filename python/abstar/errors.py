"""Error codes and the public exception type."""
from __future__ import annotations

EXIT_CODES: dict[int, str] = {
    0: "success",
    1: "fatal run error",
    2: "invalid usage or input",
    3: "strict-mode record failure",
}

ERROR_CODES: dict[str, tuple[int, str]] = {
    "INPUT_NOT_FOUND": (2, "an input file does not exist, or a string is neither a file nor a sequence"),
    "INPUT_UNREADABLE": (2, "an input file could not be read or parsed"),
    "INVALID_ARGUMENT": (2, "an argument value is invalid"),
    "OUTPUT_NAME_CONFLICT": (2, "two inputs would write the same output file"),
    "DATABASE_NOT_FOUND": (2, "the germline database does not exist"),
    "DATABASE_INVALID": (2, "the germline database is malformed"),
    "RECIPE_INVALID": (2, "the database recipe is malformed or references missing files"),
    "DATABASE_EXISTS": (2, "the output database already exists (use overwrite)"),
    "BUILTIN_NAME_CONFLICT": (2, "the database name is reserved by a built-in database"),
    "IO_ERROR": (1, "a file system read or write failed"),
    "INTERNAL": (1, "an internal error (please report it)"),
    "STRICT_RECORD_FAILURE": (3, "a record failed to annotate in strict mode"),
}


class AbstarError(Exception):
    """An abstar error with a stable machine-readable `code`."""

    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message
        self.exit_code = ERROR_CODES.get(code, (1, ""))[0]


def from_core(e: Exception) -> AbstarError:
    """Convert an `abstar._core.CoreError` (args == (code, message))."""
    code, message = e.args
    return AbstarError(code, message)
