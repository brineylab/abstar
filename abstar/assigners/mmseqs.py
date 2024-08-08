# Copyright (c) 2024 Bryan Briney
# Distributed under the terms of the MIT License.
# SPDX-License-Identifier: MIT


import logging
import os
import sys
from typing import Iterable, Optional, Union

import abutils
import polars as pl

from .assigner import AssignerBase


class MMseqs(AssignerBase):
    def __init__(
        self,
        output_directory: str,
        log_directory: str,
        germdb_name: str,
        receptor: str,
        logger: Optional[logging.Logger] = None,
        debug: bool = False,
    ) -> None:
        """
        Initialize the MMseqs assigner.
        """
        super().__init__(
            output_directory=output_directory,
            log_directory=log_directory,
            germdb_name=germdb_name,
            receptor=receptor,
            logger=logger,
            debug=debug,
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
        input_fasta, input_tsv, sequence_count = self.prepare_input_files(sequence_file)
        self.logger.info(f"found {sequence_count:,} sequences\n")
        assigned_path = self.assign_germlines(
            input_fasta=input_fasta, input_tsv=input_tsv
        )
        return assigned_path, sequence_count

    def assign_germlines(self, input_fasta: str, input_tsv: str) -> str:
        """
        V, D, J, and C germline gene segment assignment.

        Parameters
        ----------
        input_fasta : str
            The path to the input FASTA file.

        input_tsv : str
            The path to the input TSV file.

        Returns
        -------
        str
            The path to the output Parquet file.
        """
        mmseqs_format_output = "query,target,evalue,qstart,qend,qseq,nident"
        germdb_path = os.path.join(self.germdb_path, "mmseqs")
        input_df = pl.scan_csv(input_tsv, separator="\t")

        # -----------
        #   V genes
        # -----------

        self.logger.info("\n")
        self.logger.info("germline assignment:\n")
        self.logger.info("  V")
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
            sensitivity=7.5,
            # max_evalue=1.0e-6,
            format_mode=4,
            additional_cli_args="--alignment-mode 3",
            format_output=mmseqs_format_output,
            log_to=os.path.join(self.log_directory, "v_assignment.log"),
            debug=self.debug,
        )
        # read the results
        vresult_df = pl.scan_csv(  # vresult_df is a LazyFrame
            vresult_path,
            separator="\t",
            with_column_names=lambda x: [
                f"v_{_x}".replace("target", "call").replace("evalue", "support")
                for _x in x
            ],
        )
        # keep only the highest scoring assignment for each sequence
        # vresult_df = vresult_df.sort(by=["v_support"], nulls_last=True)
        vresult_df = vresult_df.sort(by=["v_nident"], descending=True, nulls_last=True)
        vresult_df = vresult_df.unique(subset=["v_query"], keep="first")

        # -----------
        #   J genes
        # -----------

        self.logger.info(" | J")
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
            additional_cli_args="--min-aln-len 12 -k 5 --alignment-mode 3",
            format_output=mmseqs_format_output,
            log_to=os.path.join(self.log_directory, "j_assignment.log"),
            debug=self.debug,
        )
        # read the results
        jresult_df = pl.scan_csv(
            jresult_path,
            separator="\t",
            with_column_names=lambda x: [
                f"j_{_x}".replace("target", "call").replace("evalue", "support")
                for _x in x
            ],
        )
        # keep only the highest scoring assignment for each sequence
        jresult_df = jresult_df.sort(by=["j_nident"], descending=True, nulls_last=True)
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

        self.logger.info(" | D")
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
        self.build_dquery_fasta(vjresult_df, dquery_path)
        # assign D genes
        abutils.tl.mmseqs_search(
            query=dquery_path,
            target=d_germdb,
            output_path=dresult_path,
            search_type=3,
            max_seqs=25,
            max_evalue=1000000.0,
            format_mode=4,
            additional_cli_args="--min-aln-len 5 -k 3 --alignment-mode 3",
            format_output=mmseqs_format_output,
            log_to=os.path.join(self.log_directory, "d_assignment.log"),
            debug=self.debug,
        )
        # read the results
        dresult_df = pl.scan_csv(
            dresult_path,
            separator="\t",
            with_column_names=lambda x: [
                f"d_{_x}".replace("target", "call").replace("evalue", "support")
                for _x in x
            ],
        )
        # keep only the highest scoring assignment for each sequence
        dresult_df = dresult_df.sort(by=["d_nident"], descending=True, nulls_last=True)
        dresult_df = dresult_df.unique(subset=["d_query"], keep="first")
        if dresult_df.select(pl.len()).collect().item() > 0:
            # join the D and VJ assignment results
            vdjresult_df = vjresult_df.join(
                dresult_df,
                left_on="v_query",
                right_on="d_query",
                how="left",
            )
        else:
            # if none of the sequences have a D gene assignment,
            # set the D gene columns to None
            vdjresult_df = vjresult_df.with_columns(
                pl.lit(None).alias("d_call"),
                pl.lit(None).alias("d_support"),
            )

        # -----------
        #   C genes
        # -----------
        c_germdb = os.path.join(germdb_path, "c")
        cquery_path = os.path.join(
            self.output_directory, f"{self.sample_name}.cquery.fasta"
        )
        # some germline databases may not have constant genes
        # so we only try to assign if the database exists
        if os.path.exists(c_germdb):
            self.logger.info(" | C")
            cresult_path = os.path.join(
                self.output_directory, f"{self.sample_name}.cresult.tsv"
            )
            if not self.debug:
                self.to_delete.extend([cresult_path, cquery_path])
            # make the input FASTA file for C gene assignment
            self.build_cquery_fasta(vjresult_df, cquery_path)
            # assign C genes
            abutils.tl.mmseqs_search(
                query=cquery_path,
                target=c_germdb,
                output_path=cresult_path,
                search_type=3,
                # max_seqs=25,
                max_evalue=10.0,
                format_mode=4,
                additional_cli_args="--min-aln-len 12 -k 5 --alignment-mode 3",
                format_output=mmseqs_format_output,
                log_to=os.path.join(self.log_directory, "c_assignment.log"),
                debug=self.debug,
            )
            # read the results
            cresult_df = pl.scan_csv(
                cresult_path,
                separator="\t",
                with_column_names=lambda x: [
                    f"c_{_x}".replace("target", "call").replace("evalue", "support")
                    for _x in x
                ],
            )
            # keep only the highest scoring assignment for each sequence
            cresult_df = cresult_df.sort(
                by=["c_nident"], descending=True, nulls_last=True
            )
            cresult_df = cresult_df.unique(subset=["c_query"], keep="first")
            if cresult_df.select(pl.len()).collect().item() > 0:
                # join the C and VDJ assignment results
                vdjcresult_df = vdjresult_df.join(
                    cresult_df,
                    left_on="v_query",
                    right_on="c_query",
                    how="left",
                )
            else:
                # if none of the sequences have a C gene assignment,
                # set the C gene columns to None
                vdjcresult_df = vdjresult_df.with_columns(
                    pl.lit(None).alias("c_call"),
                    pl.lit(None).alias("c_support"),
                )
        else:
            # if the C segment germline database does not exist,
            # set the C gene columns to None
            vdjcresult_df = vdjresult_df.with_columns(
                pl.lit(None).alias("c_call"),
                pl.lit(None).alias("c_support"),
            )

        # join the input CSV data with VDJC assignment results
        self.logger.info("\n")
        self.logger.info("combining germline assignment results\n")
        vdjcresult_df = input_df.join(
            vdjcresult_df,
            left_on="sequence_id",
            right_on="v_query",
            how="left",
        )

        # add the rev_comp column
        vdjcresult_df = vdjcresult_df.with_columns(
            (pl.col("v_qstart") > pl.col("v_qend")).alias("rev_comp")
        )

        # log "unassigned" sequences (no V or J gene assignment)
        unassigned_cols = ["sequence_id", "sequence_input"]
        unassigned = vdjcresult_df.filter(
            pl.col("v_call").is_null() | pl.col("j_call").is_null()
        )
        unassigned = unassigned.select(unassigned_cols)
        unassigned_path = os.path.join(
            self.log_directory, f"{self.sample_name}.unassigned.csv"
        )
        unassigned.sink_csv(unassigned_path)

        # write "assigned" sequence results (successful V and J gene assignment)
        assigned_cols = [
            "sequence_id",
            "sequence_input",
            "quality",
            "rev_comp",
            "v_call",
            "v_support",
            "d_call",
            "d_support",
            "j_call",
            "j_support",
            "c_call",
            "c_support",
        ]
        assigned = vdjcresult_df.filter(
            ~pl.col("v_call").is_null() & ~pl.col("j_call").is_null()
        )
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
            tab-delimited file containing the sequence ID, input sequence, and (optionally) quality score
            for each input sequence.

            .. note:
                The output format needs to be tab-separated rather than comma-separated because the comma
                is a valid character in quality scores and will obviously mess with CSV formatting.
        """
        # set up output files
        output_fasta = os.path.join(self.output_directory, f"{self.sample_name}.fasta")
        output_csv = os.path.join(self.output_directory, f"{self.sample_name}.tsv")
        if not self.debug:
            self.to_delete.extend([output_fasta, output_csv])

        # process input file
        sequence_count = 0
        with open(output_fasta, "w") as ofasta:
            with open(output_csv, "w") as ocsv:
                ocsv.write("sequence_id\tsequence_input\tquality\n")  # header
                for seq in abutils.io.parse_fastx(sequence_file):
                    qual = seq.qual if seq.qual is not None else ""
                    ofasta.write(f">{seq.id}\n{seq.sequence}\n")
                    ocsv.write(f"{seq.id}\t{seq.sequence}\t{qual}\n")
                    sequence_count += 1
        return output_fasta, output_csv, sequence_count

    def build_jquery_fasta(
        self,
        vresult_df: Union[pl.LazyFrame, pl.DataFrame],
        fasta_path: str,
    ) -> None:
        """
        Builds a FASTA file for the J gene assignment.
        The query sequence is the portion of the input sequence
        that follows the region assigned to the V gene.

        Parameters
        ----------
        vresult_df : pl.LazyFrame or pl.DataFrame
            The V gene assignment results.

        fasta_path : str
            The path to the output FASTA file.

        Returns
        -------
        None

        """
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
        self,
        vjresult_df: Union[pl.LazyFrame, pl.DataFrame],
        fasta_path: str,
    ) -> None:
        """
        Builds a FASTA file for the D gene assignment.
        The query sequence is the portion of the input sequence
        that falls between the regions assigned to V and J genes.

        Parameters
        ----------
        vjresult_df : pl.LazyFrame or pl.DataFrame
            The J gene assignment results.

        fasta_path : str
            The path to the output FASTA file.

        Returns
        -------
        None

        """
        vcols = ["v_query", "j_qstart", "j_qend", "j_qseq"]
        fastas = []

        # only assign D genes on "heavy" chains
        prefixes = ["IGH", "TRA", "TRD"]
        heavy_df = vjresult_df.filter(pl.col("v_call").str.contains_any(prefixes))
        heavy_df = heavy_df.filter(~pl.col("j_qstart").is_null())
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

    def build_cquery_fasta(
        self,
        vjresult_df: Union[pl.LazyFrame, pl.DataFrame],
        fasta_path: str,
    ) -> None:
        """
        Builds a FASTA file for the C gene assignment.
        The query sequence is the portion of the input sequence
        that follows the region assigned to the J gene.

        Parameters
        ----------
        vjresult_df : pl.LazyFrame or pl.DataFrame
            The J gene assignment results.

        fasta_path : str
            The path to the output FASTA file.

        Returns
        -------
        None

        """
        vcols = ["v_query", "j_qstart", "j_qend", "j_qseq"]
        fastas = []
        if isinstance(vjresult_df, pl.LazyFrame):
            vjresult_df = vjresult_df.collect()
        vjresult_df = vjresult_df.drop_nulls(vcols)

        # retrieve a subset of the sequence that preceeds the J alignment
        for name, start, end, seq in vjresult_df.select(
            pl.concat_list(pl.col(vcols)).alias("cq")
        )["cq"]:
            start = int(start)
            end = int(end)
            if start > end:
                query = seq[:end]
            else:
                query = seq[end:]
            fasta = f">{name}\n{query}"
            fastas.append(fasta)

        # write the FASTA file
        with open(fasta_path, "w") as f:
            f.write("\n".join(fastas))


# class MMseqs(AssignerBase):
#     def __init__(
#         self,
#         output_directory: str,
#         log_directory: str,
#         germdb_name: str,
#         receptor: str,
#         verbose: bool = False,
#         debug: bool = False,
#     ) -> None:
#         """
#         Initialize the MMseqs assigner.
#         """
#         super().__init__(
#             output_directory=output_directory,
#             log_directory=log_directory,
#             germdb_name=germdb_name,
#             receptor=receptor,
#             verbose=verbose,
#             debug=debug,
#         )

#     def __call__(self, sequence_file: str) -> str:
#         """
#         Run the MMseqs assigner.

#         Parameters
#         ----------
#         sequence_file : str
#             The path to the input file, in either FASTA or FASTQ format. Gzip-compressed files
#             are supported.

#         Returns
#         -------
#         str
#             The path to the output Parquet file.
#         """
#         self.sample_name = ".".join(
#             os.path.basename(sequence_file).rstrip(".gz").split(".")[:-1]
#         )

#         if self.verbose:
#             print("  - preparing input file")
#         input_fasta, input_tsv = self.prepare_input_files(sequence_file)
#         return self.assign_germlines(input_fasta=input_fasta, input_tsv=input_tsv)

#     def assign_germlines(self, input_fasta: str, input_tsv: str) -> str:
#         """
#         V, D and J germline gene segment assignment.

#         Parameters
#         ----------
#         input_fasta : str
#             The path to the input FASTA file.

#         input_tsv : str
#             The path to the input TSV file.

#         Returns
#         -------
#         str
#             The path to the output Parquet file.
#         """
#         mmseqs_format_output = "query,target,evalue,qstart,qend,qseq"
#         germdb_path = os.path.join(self.germdb_path, "mmseqs")
#         input_df = pl.scan_csv(input_tsv, separator="\t")

#         # -----------
#         #   V genes
#         # -----------
#         if self.verbose:
#             print("  - assigning V genes")
#         v_germdb = os.path.join(germdb_path, "v")
#         vresult_path = os.path.join(
#             self.output_directory, f"{self.sample_name}.vresult.tsv"
#         )
#         if not self.debug:
#             self.to_delete.append(vresult_path)
#         # assign V genes
#         abutils.tl.mmseqs_search(
#             query=input_fasta,
#             target=v_germdb,
#             output_path=vresult_path,
#             search_type=3,
#             max_seqs=25,
#             max_evalue=1.0e-6,
#             format_mode=4,
#             format_output=mmseqs_format_output,
#             debug=self.debug,
#         )
#         # read the results
#         vresult_df = pl.scan_csv(  # vresult_df is a LazyFrame
#             vresult_path,
#             separator="\t",
#             # truncate_ragged_lines=True,
#             with_column_names=lambda x: [
#                 f"v_{_x}".replace("target", "call") for _x in x
#             ],
#         )
#         # keep only the highest scoring assignment for each sequence
#         vresult_df = vresult_df.sort(by=["v_evalue"], nulls_last=True)
#         vresult_df = vresult_df.unique(subset=["v_query"], keep="first")

#         # -----------
#         #   J genes
#         # -----------
#         if self.verbose:
#             print("  - assigning J genes")
#         j_germdb = os.path.join(germdb_path, "j")
#         jquery_path = os.path.join(
#             self.output_directory, f"{self.sample_name}.jquery.fasta"
#         )
#         jresult_path = os.path.join(
#             self.output_directory, f"{self.sample_name}.jresult.tsv"
#         )
#         if not self.debug:
#             self.to_delete.extend([jresult_path, jquery_path])
#         # make the input FASTA file for J gene assignment
#         self.build_jquery_fasta(vresult_df, jquery_path)
#         # assign J genes
#         abutils.tl.mmseqs_search(
#             query=jquery_path,
#             target=j_germdb,
#             output_path=jresult_path,
#             search_type=3,
#             max_seqs=25,
#             max_evalue=1000.0,
#             format_mode=4,
#             additional_cli_args="--min-aln-len 12 -k 3",
#             format_output=mmseqs_format_output,
#             debug=self.debug,
#         )
#         # read the results
#         jresult_df = pl.scan_csv(
#             jresult_path,
#             separator="\t",
#             # truncate_ragged_lines=True,
#             with_column_names=lambda x: [
#                 f"j_{_x}".replace("target", "call") for _x in x
#             ],
#         )
#         # keep only the highest scoring assignment for each sequence
#         jresult_df = jresult_df.sort(by=["j_evalue"], nulls_last=True)
#         jresult_df = jresult_df.unique(subset=["j_query"], keep="first")

#         # join the V and J assignment results
#         vjresult_df = vresult_df.join(
#             jresult_df,
#             left_on="v_query",
#             right_on="j_query",
#             how="left",
#         )

#         # -----------
#         #   D genes
#         # -----------
#         if self.verbose:
#             print("  - assigning D genes")
#         d_germdb = os.path.join(germdb_path, "d")
#         dquery_path = os.path.join(
#             self.output_directory, f"{self.sample_name}.dquery.fasta"
#         )
#         dresult_path = os.path.join(
#             self.output_directory, f"{self.sample_name}.dresult.tsv"
#         )
#         if not self.debug:
#             self.to_delete.extend([dresult_path, dquery_path])
#         # make the input FASTA file for D gene assignment
#         self.build_dquery_fasta(vjresult_df, dquery_path)
#         # assign D genes
#         abutils.tl.mmseqs_search(
#             query=dquery_path,
#             target=d_germdb,
#             output_path=dresult_path,
#             search_type=3,
#             max_seqs=25,
#             max_evalue=10000.0,
#             format_mode=4,
#             additional_cli_args="--min-aln-len 5 -k 3",
#             format_output=mmseqs_format_output,
#             debug=self.debug,
#         )
#         # read the results
#         dresult_df = pl.scan_csv(
#             dresult_path,
#             separator="\t",
#             # truncate_ragged_lines=True,
#             with_column_names=lambda x: [
#                 f"d_{_x}".replace("target", "call") for _x in x
#             ],
#         )
#         # keep only the highest scoring assignment for each sequence
#         dresult_df = dresult_df.sort(by=["d_evalue"], nulls_last=True)
#         dresult_df = dresult_df.unique(subset=["d_query"], keep="first")

#         # join the D and VJ assignment results
#         vdjresult_df = vjresult_df.join(
#             dresult_df,
#             left_on="v_query",
#             right_on="d_query",
#             how="left",
#         )

#         # join the input CSV data with VDJ assignment results
#         vdjresult_df = input_df.join(
#             vdjresult_df,
#             left_on="sequence_id",
#             right_on="v_query",
#             how="left",
#         )

#         # add the is_rc column
#         vdjresult_df = vdjresult_df.with_columns(
#             (pl.col("v_qstart") > pl.col("v_qend")).alias("is_rc")
#         )

#         # log "unassigned" sequences (no V or J gene assignment)
#         unassigned_cols = ["sequence_id", "sequence_input"]
#         unassigned = vdjresult_df.filter(
#             pl.col("v_call").is_null() | pl.col("j_call").is_null()
#         )
#         unassigned = unassigned.select(unassigned_cols)
#         unassigned_path = os.path.join(
#             self.log_directory, f"{self.sample_name}.unassigned.csv"
#         )
#         unassigned.sink_csv(unassigned_path)

#         # write "assigned" sequence results (successful V and J gene assignment)
#         assigned_cols = [
#             "sequence_id",
#             "sequence_input",
#             "quality",
#             "is_rc",
#             "v_call",
#             "v_evalue",
#             "d_call",
#             "d_evalue",
#             "j_call",
#             "j_evalue",
#         ]
#         assigned = vdjresult_df.filter(
#             ~pl.col("v_call").is_null() & ~pl.col("j_call").is_null()
#         )
#         assigned = assigned.select(assigned_cols)
#         assigned_path = os.path.join(
#             self.output_directory, f"{self.sample_name}.parquet"
#         )
#         assigned.sink_parquet(assigned_path)

#         return assigned_path

#     def prepare_input_files(self, sequence_file: str) -> Iterable[str]:
#         """
#         Prepare the input file for use by the Assigner.

#         MMSeqs requires FASTA-formatted files, but the input file may be FASTQ-formatted and
#         may also be gzip-compressed. This method handles the conversion and decompression of the
#         input file.

#         Parameters
#         ----------
#         sequence_file : str
#             The path to the input file.

#         Returns
#         -------
#         Iterable[str]
#             An iterable containing the path to a FASTA-formatted input file, and the path to a 3-column
#             tab-delimited file containing the sequence ID, input sequence, and (optionally) quality score
#             for each input sequence.

#             .. note:
#                 The output format needs to be tab-separated rather than comma-separated because the comma
#                 is a valid character in quality scores and will obviously mess with CSV formatting.
#         """
#         # set up output files
#         output_fasta = os.path.join(self.output_directory, f"{self.sample_name}.fasta")
#         output_csv = os.path.join(self.output_directory, f"{self.sample_name}.tsv")
#         if not self.debug:
#             self.to_delete.extend([output_fasta, output_csv])

#         # process input file
#         with open(output_fasta, "w") as ofasta:
#             with open(output_csv, "w") as ocsv:
#                 ocsv.write("sequence_id\tsequence_input\tquality\n")  # header
#                 for seq in abutils.io.parse_fastx(sequence_file):
#                     qual = seq.qual if seq.qual is not None else ""
#                     ofasta.write(f">{seq.id}\n{seq.sequence}\n")
#                     ocsv.write(f"{seq.id}\t{seq.sequence}\t{qual}\n")
#         return output_fasta, output_csv

#     def build_jquery_fasta(
#         self,
#         vresult_df: Union[pl.LazyFrame, pl.DataFrame],
#         fasta_path: str,
#     ) -> None:
#         """ """
#         vcols = ["v_query", "v_qstart", "v_qend", "v_qseq"]
#         fastas = []
#         if isinstance(vresult_df, pl.LazyFrame):
#             vresult_df = vresult_df.collect()

#         # retrieve a subset of the sequence that follows the V alignment
#         for name, start, end, seq in vresult_df.select(
#             pl.concat_list(pl.col(vcols)).alias("jq")
#         )["jq"]:
#             start = int(start)
#             end = int(end)
#             if start > end:
#                 fasta = f">{name}\n{seq[:end]}"
#             else:
#                 fasta = f">{name}\n{seq[end:]}"
#             fastas.append(fasta)

#         # write the FASTA file
#         with open(fasta_path, "w") as f:
#             f.write("\n".join(fastas))

#     def build_dquery_fasta(
#         self,
#         vjresult_df: Union[pl.LazyFrame, pl.DataFrame],
#         fasta_path: str,
#     ) -> None:
#         """ """
#         vcols = ["v_query", "j_qstart", "j_qend", "j_qseq"]
#         fastas = []

#         # only assign D genes on "heavy" chains
#         prefixes = ["IGH", "TRA", "TRD"]
#         heavy_df = vjresult_df.filter(pl.col("v_call").str.contains_any(prefixes))
#         heavy_df = heavy_df.filter(~pl.col("j_qstart").is_null())
#         if isinstance(heavy_df, pl.LazyFrame):
#             heavy_df = heavy_df.collect()

#         # retrieve a subset of the sequence that preceeds the J alignment
#         for name, start, end, seq in heavy_df.select(
#             pl.concat_list(pl.col(vcols)).alias("dq")
#         )["dq"]:
#             start = int(start)
#             end = int(end)
#             if start > end:
#                 fasta = f">{name}\n{seq[start:]}"
#             else:
#                 fasta = f">{name}\n{seq[:start]}"
#             fastas.append(fasta)

#         # write the FASTA file
#         with open(fasta_path, "w") as f:
#             f.write("\n".join(fastas))
