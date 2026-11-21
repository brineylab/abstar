//! Shared test helpers: the built-in human database and generated reads with their truth
//! (generator dev-dependency; the library never uses it).
#![allow(dead_code)]
use abstar_annotate::align::Scoring;
use abstar_annotate::assign::{assign_c, assign_j, call_v, resolve_junction, CHit, JHit, Junction, Oriented, SegmentIndex, VHit};
use abstar_generate::{bind, generate_in_memory, load_model, preset, GenerateParams};
use abstar_germline::Database;
use abstar_germline::{Allele, Segment};
use arrow_array::{Array, BooleanArray, Float64Array, Int64Array, ListArray, RecordBatch, StringArray};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

pub fn package_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../python/abstar")
}

/// The built-in `human` germline database (IG and TR).
pub fn human_db() -> Database {
    abstar_germline::load(&package_root().join("germline_dbs/human")).unwrap()
}

/// A database's alleles with their V, D, J and C indexes.
pub struct Ref {
    pub alleles: Vec<Allele>,
    pub v: SegmentIndex,
    pub d: SegmentIndex,
    pub j: SegmentIndex,
    pub c: SegmentIndex,
}

/// Every call for one read: V, J, the junction and C.
pub struct Calls {
    pub read: Oriented,
    pub v: VHit,
    pub j: JHit,
    pub junction: Junction,
    pub c: Option<CHit>,
}

impl Ref {
    pub fn new(alleles: Vec<Allele>) -> Ref {
        let v = SegmentIndex::new(&alleles, Segment::V);
        let d = SegmentIndex::new(&alleles, Segment::D);
        let j = SegmentIndex::new(&alleles, Segment::J);
        let c = SegmentIndex::new(&alleles, Segment::C);
        Ref { alleles, v, d, j, c }
    }
    pub fn call(&self, read: &[u8]) -> Option<(Oriented, VHit)> {
        call_v(read, &self.alleles, &self.v, &self.j, &Scoring::default())
    }
    /// The V call, then J assignment after it.
    pub fn call_j(&self, read: &[u8]) -> Option<JHit> {
        let (o, v) = self.call(read)?;
        assign_j(&o, &v, &self.alleles, &self.j, &Scoring::default())
    }
    /// V, J, then the junction (a D searched when the V locus has D alleles) and C.
    pub fn call_all(&self, read: &[u8]) -> Option<Calls> {
        let (o, v) = self.call(read)?;
        let j = assign_j(&o, &v, &self.alleles, &self.j, &Scoring::default())?;
        let has_d = self.d.loci.contains(&self.alleles[v.primary].locus);
        let junction = resolve_junction(&o, &v, &j, &self.alleles, &self.d, has_d);
        let c = assign_c(&o, &j, &self.alleles, &self.c);
        Some(Calls { read: o, v, j, junction, c })
    }
    pub fn idx(&self, name: &str) -> usize {
        self.alleles.iter().position(|a| a.name == name).unwrap_or_else(|| panic!("no allele {name}"))
    }
    pub fn names(&self, ids: &[usize]) -> BTreeSet<String> {
        ids.iter().map(|&i| self.alleles[i].name.clone()).collect()
    }
}

/// The built-in `human` database, indexed once.
pub fn human() -> &'static Ref {
    static R: OnceLock<Ref> = OnceLock::new();
    R.get_or_init(|| Ref::new(human_db().alleles))
}

/// Generated reads and their truth table (one Arrow batch, one row per read).
pub struct Generated(pub RecordBatch);

impl Generated {
    pub fn len(&self) -> usize {
        self.0.num_rows()
    }
    fn col<T: 'static>(&self, name: &str) -> &T {
        let c = self.0.column_by_name(name).unwrap_or_else(|| panic!("no truth column {name}"));
        c.as_any().downcast_ref::<T>().unwrap_or_else(|| panic!("truth column {name} has another type"))
    }
    pub fn str(&self, name: &str, i: usize) -> Option<&str> {
        let a = self.col::<StringArray>(name);
        a.is_valid(i).then(|| a.value(i))
    }
    pub fn int(&self, name: &str, i: usize) -> Option<i64> {
        let a = self.col::<Int64Array>(name);
        a.is_valid(i).then(|| a.value(i))
    }
    pub fn bool(&self, name: &str, i: usize) -> Option<bool> {
        let a = self.col::<BooleanArray>(name);
        a.is_valid(i).then(|| a.value(i))
    }
    pub fn list(&self, name: &str, i: usize) -> Vec<String> {
        let a = self.col::<ListArray>(name);
        let v = a.value(i);
        let s = v.as_any().downcast_ref::<StringArray>().unwrap();
        (0..s.len()).map(|j| s.value(j).to_string()).collect()
    }
    pub fn float(&self, name: &str, i: usize) -> Option<f64> {
        let a = self.col::<Float64Array>(name);
        a.is_valid(i).then(|| a.value(i))
    }
    /// A list<string> cell; `None` when null.
    pub fn opt_list(&self, name: &str, i: usize) -> Option<Vec<String>> {
        self.col::<ListArray>(name).is_valid(i).then(|| self.list(name, i))
    }
    /// A list<int64> cell; empty when null.
    pub fn ints(&self, name: &str, i: usize) -> Vec<i64> {
        let a = self.col::<ListArray>(name);
        if !a.is_valid(i) {
            return Vec::new();
        }
        a.value(i).as_any().downcast_ref::<Int64Array>().unwrap().values().to_vec()
    }
    pub fn sequence(&self, i: usize) -> &str {
        self.str("sequence", i).unwrap()
    }
}

/// `n` reads from the built-in `human-igh` model on the built-in `human` database, with the
/// `preset`'s parameters, `seed`, then `edit` applied.
pub fn generated(seed: u64, preset_name: &str, n: usize, edit: impl FnOnce(&mut GenerateParams)) -> Generated {
    generated_inverting(seed, preset_name, n, 0.0, edit)
}

/// `generated` with the model's `d_inversion` set to `d_inversion`.
pub fn generated_inverting(seed: u64, preset_name: &str, n: usize, d_inversion: f64, edit: impl FnOnce(&mut GenerateParams)) -> Generated {
    let root = package_root();
    let mut model = load_model(&root.join("models/human-igh/model.toml")).unwrap();
    model.d_inversion = d_inversion;
    let bound = bind(model, &human_db()).unwrap();
    let mut params = preset(preset_name).unwrap();
    params.n = n;
    params.seed = seed;
    params.batch_size = n.max(1);
    edit(&mut params);
    let mut batches = generate_in_memory(&bound, &params).unwrap();
    assert_eq!(batches.len(), 1);
    Generated(batches.pop().unwrap())
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            o => o,
        })
        .collect()
}

/// Deterministic pseudo-random ACGT (xorshift64).
pub fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
    let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x >> 33) as usize % 4]
        })
        .collect()
}
