#!/usr/bin/env python
# filename: mmseqs.py

#
# Copyright (c) 2024 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


import os
from typing import Iterable, Union

import abutils
import polars as pl

from .assigner import AssignerBase


class MMseqs(AssignerBase):
    def __init__(self, output_directory: str, germdb_name: str, receptor: str) -> None:
        """
        Initialize the MMseqs assigner.
        """
        super().__init__(
            output_directory=output_directory,
            germdb_name=germdb_name,
            receptor=receptor,
        )

    def __call__(self, sequence_file: str) -> str:
        """
        Run the MMseqs assigner.

        Parameters
        ----------
        sequence_file : str
            The path to the input file, in either FASTA or FASTQ format. Gzip-compressed files
            are supported.

        Returns
        -------
        str
            The path to the output Parquet file.
        """
        self.sample_name = ".".join(
            os.path.basename(sequence_file).rstrip(".gz").split(".")[:-1]
        )
        input_fasta, input_csv = self.prepare_input_files(sequence_file)
        return self.assign_germlines(input_fasta, input_csv)

    def assign_germlines(self, input_fasta: str, input_csv: str) -> str:
        """
        V, D and J germline gene segment assignment.

        Parameters
        ----------
        input_fasta : str
            The path to the input FASTA file.

        input_csv : str
            The path to the input CSV file.

        Returns
        -------
        str
            The path to the output Parquet file.
        """
        mmseqs_format_output = "query,target,evalue,qstart,qend,qseq"
        germdb_path = os.path.join(self.germdb_path, "mmseqs")
        input_df = pl.scan_csv(input_csv)

        # -----------
        #   V genes
        # -----------
        v_germdb = os.path.join(germdb_path, "v")
        vresult_path = os.path.join(
            self.output_directory, f"{self.sample_name}.vresult.tsv"
        )
        if not self.debug:
            self.to_delete.append(vresult_path)
        # assign V genes
        abutils.tl.mmseqs_search(
            query=input_fasta,
            target=v_germdb,
            output_path=vresult_path,
            search_type=3,
            max_seqs=25,
            max_evalue=1.0e-6,
            format_mode=4,
            format_output=mmseqs_format_output,
            debug=self.debug,
        )
        # read the results
        vresult_df = pl.scan_csv(  # vresult_df is a LazyFrame
            vresult_path,
            separator="\t",
            with_column_names=lambda x: [
                f"v_{_x}".replace("target", "call") for _x in x
            ],
        )
        # keep only the highest scoring assignment for each sequence
        vresult_df = vresult_df.sort(by=["v_evalue"], nulls_last=True)
        vresult_df = vresult_df.unique(subset=["v_query"], keep="first")

        # -----------
        #   J genes
        # -----------
        j_germdb = os.path.join(germdb_path, "j")
        jquery_path = os.path.join(
            self.output_directory, f"{self.sample_name}.jquery.fasta"
        )
        jresult_path = os.path.join(
            self.output_directory, f"{self.sample_name}.jresult.tsv"
        )
        if not self.debug:
            self.to_delete.extend([jresult_path, jquery_path])
        # make the input FASTA file for J gene assignment
        self.build_jquery_fasta(vresult_df, jquery_path)
        # assign J genes
        abutils.tl.mmseqs_search(
            query=jquery_path,
            target=j_germdb,
            output_path=jresult_path,
            search_type=3,
            max_seqs=25,
            max_evalue=1000.0,
            format_mode=4,
            additional_cli_args="--min-aln-len 16 -k 6",
            format_output=mmseqs_format_output,
            debug=self.debug,
        )
        # read the results
        jresult_df = pl.scan_csv(
            jresult_path,
            separator="\t",
            with_column_names=lambda x: [
                f"j_{_x}".replace("target", "call") for _x in x
            ],
        )
        # keep only the highest scoring assignment for each sequence
        jresult_df = jresult_df.sort(by=["j_evalue"], nulls_last=True)
        jresult_df = jresult_df.unique(subset=["j_query"], keep="first")

        # join the V and J assignment results
        vjresult_df = vresult_df.join(
            jresult_df,
            left_on="v_query",
            right_on="j_query",
            how="left",
        )

        # -----------
        #   D genes
        # -----------
        d_germdb = os.path.join(germdb_path, "d")
        dquery_path = os.path.join(
            self.output_directory, f"{self.sample_name}.dquery.fasta"
        )
        dresult_path = os.path.join(
            self.output_directory, f"{self.sample_name}.dresult.tsv"
        )
        if not self.debug:
            self.to_delete.extend([dresult_path, dquery_path])
        # make the input FASTA file for D gene assignment
        self.build_dquery_fasta(jresult_df, dquery_path)
        # assign D genes
        abutils.tl.mmseqs_search(
            query=dquery_path,
            target=d_germdb,
            output_path=dresult_path,
            search_type=3,
            max_seqs=25,
            max_evalue=10000.0,
            format_mode=4,
            additional_cli_args="--min-aln-len 5 -k 3",
            format_output=mmseqs_format_output,
            debug=self.debug,
        )
        # read the results
        dresult_df = pl.scan_csv(
            dresult_path,
            separator="\t",
            with_column_names=lambda x: [
                f"d_{_x}".replace("target", "call") for _x in x
            ],
        )
        # keep only the highest scoring assignment for each sequence
        dresult_df = dresult_df.sort(by=["d_evalue"], nulls_last=True)
        dresult_df = dresult_df.unique(subset=["d_query"], keep="first")

        # join the D and VJ assignment results
        vdjresult_df = vjresult_df.join(
            dresult_df,
            left_on="v_query",
            right_on="d_query",
            how="left",
        )

        # join the input CSV data with VDJ assignment results
        vdjresult_df = input_df.join(
            vdjresult_df,
            left_on="sequence_id",
            right_on="v_query",
            how="left",
        )

        # log "unassigned" sequences (no V gene assignment)
        unassigned_cols = ["sequence_id", "sequence_input"]
        unassigned = vdjresult_df.filter(pl.col("v_call").is_null())
        unassigned = unassigned.select(unassigned_cols)
        unassigned_path = os.path.join(
            self.log_directory, f"{self.sample_name}.unassigned.csv"
        )
        unassigned.sink_csv(unassigned_path)

        # write "assigned" sequence results (successful V gene assignment)
        assigned_cols = [
            "sequence_id",
            "sequence_input",
            "quality",
            "is_rc",
            "v_call",
            "v_evalue",
            "d_call",
            "d_evalue",
            "j_call",
            "j_evalue",
        ]
        assigned = vdjresult_df.filter(~pl.col("v_call").is_null())
        assigned = assigned.select(assigned_cols)
        assigned_path = os.path.join(
            self.output_directory, f"{self.sample_name}.parquet"
        )
        assigned.sink_parquet(assigned_path)

        return assigned_path

    def prepare_input_files(self, sequence_file: str) -> Iterable[str]:
        """
        Prepare the input file for use by the Assigner.

        MMSeqs requires FASTA-formatted files, but the input file may be FASTQ-formatted and
        may also be gzip-compressed. This method handles the conversion and decompression of the
        input file.

        Parameters
        ----------
        sequence_file : str
            The path to the input file.

        Returns
        -------
        Iterable[str]
            An iterable containing the path to a FASTA-formatted input file, and the path to a 3-column
            CSV file containing the sequence ID, input sequence, and (optionally) quality score
            for each input sequence.
        """
        # set up output files
        output_fasta = os.path.join(self.output_directory, f"{self.sample_name}.fasta")
        output_csv = os.path.join(self.output_directory, f"{self.sample_name}.csv")
        if not self.debug:
            self.to_delete.extend([output_fasta, output_csv])

        # process input file
        with open(output_fasta, "w") as ofasta:
            with open(output_csv, "w") as ocsv:
                ocsv.write("sequence_id,sequence_input,quality\n")  # header
                for seq in abutils.io.parse_fastx(sequence_file):
                    qual = seq.qual if seq.qual is not None else ""
                    ofasta.write(f">{seq.id}\n{seq.sequence}\n")
                    ocsv.write(f"{seq.id},{seq.sequence},{qual}\n")
        return output_fasta, output_csv

    def build_jquery_fasta(
        self,
        vresult_df: Union[pl.LazyFrame, pl.DataFrame],
        fasta_path: str,
    ) -> None:
        """ """
        vcols = ["v_query", "v_qstart", "v_qend", "v_qseq"]
        fastas = []
        if isinstance(vresult_df, pl.LazyFrame):
            vresult_df = vresult_df.collect()

        # retrieve a subset of the sequence that follows the V alignment
        for name, start, end, seq in vresult_df.select(
            pl.concat_list(pl.col(vcols)).alias("jq")
        )["jq"]:
            start = int(start)
            end = int(end)
            if start > end:
                fasta = f">{name}\n{seq[:end]}"
            else:
                fasta = f">{name}\n{seq[end:]}"
            fastas.append(fasta)

        # write the FASTA file
        with open(fasta_path, "w") as f:
            f.write("\n".join(fastas))

    def build_dquery_fasta(
        vjresult_df: Union[pl.LazyFrame, pl.DataFrame],
        fasta_path: str,
    ) -> None:
        """ """
        vcols = ["v_query", "j_qstart", "j_qend", "j_qseq"]
        fastas = []

        # only assign D genes on "heavy" chains
        prefixes = ["IGH", "TRA", "TRD"]
        heavy_df = vjresult_df.filter(pl.col("v_call").str.contains_any(prefixes))
        if isinstance(heavy_df, pl.LazyFrame):
            heavy_df = heavy_df.collect()

        # retrieve a subset of the sequence that preceeds the J alignment
        for name, start, end, seq in heavy_df.select(
            pl.concat_list(pl.col(vcols)).alias("dq")
        )["dq"]:
            start = int(start)
            end = int(end)
            if start > end:
                fasta = f">{name}\n{seq[start:]}"
            else:
                fasta = f">{name}\n{seq[:start]}"
            fastas.append(fasta)

        # write the FASTA file
        with open(fasta_path, "w") as f:
            f.write("\n".join(fastas))
