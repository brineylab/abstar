//! Generates a record struct, `arrow_schema()` and `to_record_batch()` from a fields TOML file.
//! Meant for `build.rs`; the generated code expects `arrow_array`/`arrow_schema` items in scope
//! (`ArrayRef`, `RecordBatch`, `DataType`, `Field`, `Schema`, `SchemaRef`, the builders, `HashMap`, `Arc`).
use std::collections::HashSet;
use std::fmt::Write;
use std::path::Path;

/// Reads `toml_path`, writes generated Rust to `out_path`, and emits cargo rerun directives.
/// Panics (failing the build) on invalid input.
pub fn generate(toml_path: &Path, out_path: &Path, record_type: &str, const_name: &str) {
    println!("cargo:rerun-if-changed={}", toml_path.display());
    let doc: toml::Table = std::fs::read_to_string(toml_path)
        .unwrap_or_else(|e| panic!("read {}: {e}", toml_path.display()))
        .parse()
        .expect("fields toml is not valid TOML");
    let fields = doc["field"].as_array().expect("fields toml needs [[field]] entries");

    let (mut st, mut sch, mut builders, mut appends, mut cols) = <(String, String, String, String, String)>::default();
    let mut seen = HashSet::new();
    for f in fields {
        let s = |k: &str| f.get(k).and_then(|v| v.as_str()).unwrap_or_else(|| panic!("field missing `{k}`: {f}"));
        let (name, ty) = (s("name"), s("type"));
        let nullable = f["nullable"].as_bool().expect("`nullable` must be a bool");
        assert!(seen.insert(name.to_string()), "duplicate field name `{name}`");
        assert!(
            matches!(s("origin"), "airr-required" | "airr-optional" | "abstar"),
            "field `{name}`: unknown origin"
        );
        if s("origin") == "abstar" {
            let d = s("decision");
            assert!(
                d.len() == 4 && d.bytes().all(|b| b.is_ascii_digit()),
                "abstar field `{name}` needs a 4-digit `decision`, got `{d}`"
            );
        }
        // (rust type, arrow type, builder, append expression)
        let (rt, at, b, ap) = match ty {
            "string" => ("String", "Utf8", "StringBuilder::new()", "append_value(r.NAME.as_str())"),
            "bool" => ("bool", "Boolean", "BooleanBuilder::new()", "append_option(r.NAME)"),
            "int64" => ("i64", "Int64", "Int64Builder::new()", "append_option(r.NAME)"),
            "float64" => ("f64", "Float64", "Float64Builder::new()", "append_option(r.NAME)"),
            "list<string>" => ("Vec<String>", "LIST", "ListBuilder::new(StringBuilder::new())", "LIST"),
            "list<int64>" => ("Vec<i64>", "LISTI", "ListBuilder::new(Int64Builder::new())", "LIST"),
            _ => panic!("field `{name}`: unknown type `{ty}`"),
        };
        let dt = if at == "LIST" {
            "DataType::List(Arc::new(Field::new(\"item\", DataType::Utf8, true)))".to_string()
        } else if at == "LISTI" {
            "DataType::List(Arc::new(Field::new(\"item\", DataType::Int64, true)))".to_string()
        } else {
            format!("DataType::{at}")
        };
        let rust_ty = match (ty, nullable) {
            ("string", false) => rt.to_string(),
            _ => format!("Option<{rt}>"),
        };
        let ap = match (ty, nullable) {
            ("string", true) => "append_option(r.NAME.as_deref())".to_string(),
            ("list<string>" | "list<int64>", _) => "LISTAPPEND".to_string(),
            _ => ap.replace("NAME", name),
        }
        .replace("NAME", name);
        let ap = if ap == "LISTAPPEND" {
            format!(
                "match &r.{name} {{ Some(v) => {{ for x in v {{ b.values().append_value(x); }} b.append(true) }} None => b.append(false) }}"
            )
        } else {
            format!("b.{ap}")
        };
        let desc = s("description");
        writeln!(st, "    /// {desc}\n    pub {name}: {rust_ty},").unwrap();
        writeln!(sch, "        Field::new(\"{name}\", {dt}, {nullable}),").unwrap();
        writeln!(builders, "    let mut {name} = {b};").unwrap();
        writeln!(appends, "        {{ let b = &mut {name}; {ap}; }}").unwrap();
        writeln!(cols, "        Arc::new({name}.finish()) as ArrayRef,").unwrap();
    }

    let airr_version = doc["airr_version"].as_str().expect("airr_version");
    let abs = std::fs::canonicalize(toml_path).expect("canonicalize toml path");
    let code = format!(
        r#"pub const {const_name}: &str = include_str!({abs:?});

#[derive(Default, Clone, Debug, PartialEq)]
pub struct {record_type} {{
{st}}}

pub fn arrow_schema() -> SchemaRef {{
    let md: HashMap<String, String> = [
        ("airr_version".to_string(), "{airr_version}".to_string()),
        ("abstar_version".to_string(), env!("CARGO_PKG_VERSION").to_string()),
    ]
    .into();
    Arc::new(Schema::new_with_metadata(vec![
{sch}    ], md))
}}

pub fn to_record_batch(records: &[{record_type}]) -> RecordBatch {{
{builders}    for r in records {{
{appends}    }}
    RecordBatch::try_new(arrow_schema(), vec![
{cols}    ])
    .expect("generated columns match generated schema")
}}
"#
    );
    std::fs::write(out_path, code).unwrap();
}
