//! FASTA/FASTQ reading (gzip detected by magic bytes) and output-stem naming.
//! Leading BOM and blank lines are skipped; spaces and tabs are stripped from sequences.
use crate::error::RunError;
use flate2::bufread::MultiGzDecoder;
use needletail::errors::ParseErrorKind;
use needletail::FastxReader;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;

pub enum Input {
    File(PathBuf),
    Stdin,
    Records(Vec<(String, String)>),
}

impl Input {
    /// `source_file` value / run.json source: path as given, `-` for stdin, `records` otherwise.
    pub fn source(&self) -> String {
        match self {
            Input::File(p) => p.to_string_lossy().into_owned(),
            Input::Stdin => "-".into(),
            Input::Records(_) => "records".into(),
        }
    }

    /// Output file stem: file name minus `.gz`, then minus one FASTA/FASTQ extension.
    pub fn stem(&self) -> String {
        match self {
            Input::Stdin => "stdin".into(),
            Input::Records(_) => "records".into(),
            Input::File(p) => {
                let name = p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                let n = name.strip_suffix(".gz").unwrap_or(&name);
                let s = [".fasta", ".fa", ".fna", ".fastq", ".fq"]
                    .iter()
                    .find_map(|e| n.strip_suffix(e))
                    .unwrap_or(n);
                if s.is_empty() { name } else { s.into() }
            }
        }
    }
}

/// `(sequence_id, sequence, valid_utf8)`: a record whose ID or sequence bytes are not
/// UTF-8 carries lossily decoded strings and `false` (it becomes `invalid_encoding`).
pub(crate) type RawRecord = (String, String, bool);

fn decode(b: Vec<u8>) -> (String, bool) {
    match String::from_utf8(b) {
        Ok(s) => (s, true),
        Err(e) => (String::from_utf8_lossy(e.as_bytes()).into_owned(), false),
    }
}

/// Consumes a leading UTF-8 BOM and leading whitespace (blank lines, CRLF or not).
fn skip_preamble(r: &mut impl BufRead) -> std::io::Result<()> {
    let mut first = true;
    loop {
        let buf = r.fill_buf()?;
        let bom = if first && buf.starts_with(b"\xef\xbb\xbf") { 3 } else { 0 };
        first = false;
        let n = bom + buf[bom..].iter().take_while(|b| b.is_ascii_whitespace()).count();
        let done = n < buf.len() || buf.is_empty();
        r.consume(n);
        if done {
            return Ok(());
        }
    }
}

pub(crate) enum Source {
    Fastx { reader: Option<Box<dyn FastxReader>>, label: String, n: u64 },
    Mem(std::vec::IntoIter<(String, String)>),
}

impl Source {
    pub fn open(input: &Input) -> Result<Self, RunError> {
        let label = input.source();
        let raw: Box<dyn Read + Send> = match input {
            Input::Records(v) => return Ok(Source::Mem(v.clone().into_iter())),
            Input::Stdin => Box::new(std::io::stdin()),
            Input::File(p) => Box::new(std::fs::File::open(p).map_err(|e| {
                if e.kind() == std::io::ErrorKind::NotFound {
                    RunError::InputNotFound(label.clone())
                } else {
                    RunError::InputUnreadable(format!("{label}: {e}"))
                }
            })?),
        };
        let unreadable = |e: &dyn std::fmt::Display| RunError::InputUnreadable(format!("{label}: {e}"));
        if [".bz2", ".xz", ".zst"].iter().any(|x| label.to_ascii_lowercase().ends_with(x)) {
            return Err(unreadable(&"compression not supported: use gzip"));
        }
        let mut br = BufReader::new(raw);
        let gz = br.fill_buf().map_err(|e| unreadable(&e))?.starts_with(&[0x1f, 0x8b]);
        let mut rd: Box<dyn BufRead + Send> =
            if gz { Box::new(BufReader::new(MultiGzDecoder::new(br))) } else { Box::new(br) };
        skip_preamble(&mut rd).map_err(|e| unreadable(&e))?;
        let reader = match needletail::parse_fastx_reader(rd) {
            Ok(r) => Some(r),
            Err(e) if e.kind == ParseErrorKind::EmptyFile => None,
            Err(e) => return Err(unreadable(&e)),
        };
        Ok(Source::Fastx { reader, label, n: 0 })
    }

    /// Up to `max` records; empty means end of input.
    pub fn next_batch(&mut self, max: usize) -> Result<Vec<RawRecord>, RunError> {
        match self {
            Source::Mem(it) => Ok(it.take(max).map(|(id, seq)| (id, seq, true)).collect()),
            Source::Fastx { reader: None, .. } => Ok(vec![]),
            Source::Fastx { reader: Some(r), label, n } => {
                let mut out = Vec::with_capacity(max.min(65536));
                while out.len() < max {
                    let Some(rec) = r.next() else { break };
                    *n += 1;
                    let rec = rec.map_err(|e| RunError::InputUnreadable(format!("{label}: record {n}: {e}")))?;
                    let h = rec.id();
                    let h = &h[..h.iter().position(u8::is_ascii_whitespace).unwrap_or(h.len())];
                    let mut seq = rec.seq().into_owned();
                    seq.retain(|b| *b != b' ' && *b != b'\t'); // whitespace is not sequence content
                    let ((id, ok_id), (seq, ok_seq)) = (decode(h.to_vec()), decode(seq));
                    out.push((id, seq, ok_id && ok_seq));
                }
                Ok(out)
            }
        }
    }
}
