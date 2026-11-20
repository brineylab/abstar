//! Generates `Record`, `arrow_schema()` and `to_record_batch()` from `fields.toml`.
use std::collections::HashSet;
use std::fmt::Write;

fn main() {
    println!("cargo:rerun-if-changed=fields.toml");
    println!("cargo:rerun-if-changed=build.rs");
    let doc: toml::Table = std::fs::read_to_string("fields.toml")
        .expect("read fields.toml")
        .parse()
        .expect("fields.toml is not valid TOML");
    let fields = doc["field"].as_array().expect("fields.toml needs [[field]] entries");

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
            _ => panic!("field `{name}`: unknown type `{ty}`"),
        };
        let dt = if at == "LIST" {
            "DataType::List(Arc::new(Field::new(\"item\", DataType::Utf8, true)))".to_string()
        } else {
            format!("DataType::{at}")
        };
        let rust_ty = match (ty, nullable) {
            ("string", false) => rt.to_string(),
            _ => format!("Option<{rt}>"),
        };
        let ap = match (ty, nullable) {
            ("string", true) => "append_option(r.NAME.as_deref())".to_string(),
            ("list<string>", _) => "LISTAPPEND".to_string(),
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
        let doc = s("description");
        writeln!(st, "    /// {doc}\n    pub {name}: {rust_ty},").unwrap();
        writeln!(sch, "        Field::new(\"{name}\", {dt}, {nullable}),").unwrap();
        writeln!(builders, "    let mut {name} = {b};").unwrap();
        writeln!(appends, "        {{ let b = &mut {name}; {ap}; }}").unwrap();
        writeln!(cols, "        Arc::new({name}.finish()) as ArrayRef,").unwrap();
    }

    let airr_version = doc["airr_version"].as_str().expect("airr_version");
    let code = format!(
        r#"pub const FIELDS_TOML: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/fields.toml"));

#[derive(Default, Clone, Debug, PartialEq)]
pub struct Record {{
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

pub fn to_record_batch(records: &[Record]) -> RecordBatch {{
{builders}    for r in records {{
{appends}    }}
    RecordBatch::try_new(arrow_schema(), vec![
{cols}    ])
    .expect("generated columns match generated schema")
}}
"#
    );
    let out = std::path::Path::new(&std::env::var("OUT_DIR").unwrap()).join("schema_generated.rs");
    std::fs::write(out, code).unwrap();
}
