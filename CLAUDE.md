# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

abstar is a VDJ assignment and antibody/TCR sequence annotation tool. It performs germline gene assignment using MMseqs2 and detailed sequence annotation including mutations, indels, regions (CDR/FWR), and productivity assessment. Scalable from single sequences to billions.

## Common Commands

### Running Tests
```bash
# Run full test suite
pytest

# Run a single test file
pytest abstar/tests/test_regions.py

# Run a specific test
pytest abstar/tests/test_regions.py::test_get_region_sequence_fwr1

# Run with verbose output
pytest -v
```

### Installation
```bash
pip install -e .  # Development install
pip install abstar  # Production install
```

### CLI Usage
```bash
# Annotate sequences
abstar run path/to/sequences.fasta path/to/project_directory

# Build custom germline database
abstar build_germline_database <name> -f <fasta_files> -j <json_files>
```

### Python API
```python
import abstar

# Annotate sequences (returns Sequence objects when no project_path)
sequences = abstar.run("sequences.fasta")

# Annotate with output files
abstar.run("sequences.fasta", "project_dir/", output_format=["airr", "parquet"])

# Return as polars DataFrame
df = abstar.run("sequences.fasta", as_dataframe=True)
```

## Architecture

### Core Pipeline Flow
1. **Input Processing** (`core/abstar.py`): Accepts FASTA/FASTQ files, directories, or `Sequence` objects
2. **Preprocessing** (`preprocess/merging.py`): Optional paired-end read merging via fastp
3. **VDJC Assignment** (`assigners/mmseqs.py`): MMseqs2-based germline gene assignment for V, D, J, and C segments
4. **Annotation** (`annotation/annotator.py`): Parallel annotation of assigned sequences

### Key Modules

**`core/abstar.py`** - Main `run()` function orchestrating the pipeline. Handles input processing, directory setup, and parallel annotation via ProcessPoolExecutor.

**`assigners/mmseqs.py`** - `MMseqs` class performs sequential V→J→D→C gene assignment. Each step uses the previous alignment to determine query regions for the next segment.

**`annotation/annotator.py`** - `annotate_single_sequence()` performs detailed annotation:
- Germline realignment with semiglobal/local alignment
- Junction/CDR3 identification using FR3/FR4 boundaries
- Region extraction (FWR1-4, CDR1-3)
- Mutation and indel annotation
- Productivity assessment

**`annotation/antibody.py`** - `Antibody` dataclass stores all annotation fields. Includes `LoggingMixin` for debug logging.

### Module Namespaces
- `abstar.gl` - Germline database functions
- `abstar.pp` - Preprocessing (read merging)
- `abstar.tl` - Tools (UMI parsing, germline database building)

### Data Flow
```
Input → MMseqs Assignment → Parquet (temp) → Parallel Annotation → AIRR TSV / Parquet output
```

### Dependencies
- **abutils**: Provides `Sequence` objects, alignment functions, I/O utilities
- **polars**: DataFrame operations (LazyFrame for memory efficiency)
- **parasail**: Sequence alignment
- **MMseqs2**: External tool for fast sequence search (must be installed separately)

### Germline Databases
Built-in databases: `human`, `mouse`, `macaque`, `humouse` (mixed species). Custom databases can be built via `abstar.tl.build_germline_database()`.

### Output Schema
AIRR-compliant output defined in `annotation/schema.py`. Key fields include sequence_id, v/d/j/c_call, junction_aa, cdr3, mutations, productivity status.

## Testing Conventions
- Tests located in `abstar/tests/`
- Use pytest fixtures for shared test data
- Test files follow `test_<module>.py` naming
