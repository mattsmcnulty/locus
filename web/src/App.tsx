import { useEffect, useState } from "react";
import {
  api,
  classifyQuery,
  type Overview,
  type PgxResult,
  type SqlResult,
  type Variant,
  type VariantPage,
} from "./api";

type Tab = "search" | "clinical" | "pgx" | "sql";

export function App() {
  const [overview, setOverview] = useState<Overview | null>(null);
  const [err, setErr] = useState<string | null>(null);
  const [tab, setTab] = useState<Tab>("search");

  useEffect(() => {
    api.overview().then(setOverview).catch((e) => setErr(String(e.message ?? e)));
  }, []);

  return (
    <div className="app">
      <header>
        <h1>
          <span className="logo">◈</span> Locus
        </h1>
        <p className="tagline">Explore your genome locally</p>
      </header>

      {err && <div className="banner error">{err}</div>}
      {overview && <OverviewBar o={overview} />}

      <nav className="tabs">
        {(["search", "clinical", "pgx", "sql"] as Tab[]).map((t) => (
          <button key={t} className={tab === t ? "active" : ""} onClick={() => setTab(t)}>
            {t === "pgx" ? "Pharmacogenomics" : t[0].toUpperCase() + t.slice(1)}
          </button>
        ))}
      </nav>

      <main>
        {tab === "search" && <SearchView />}
        {tab === "clinical" && <ClinicalView />}
        {tab === "pgx" && <PgxView />}
        {tab === "sql" && <SqlView />}
      </main>

      <footer>
        Not medical advice — for personal exploration. Confirm health-relevant findings with a clinician.
      </footer>
    </div>
  );
}

function OverviewBar({ o }: { o: Overview }) {
  const stat = (label: string, value: number | string) => (
    <div className="stat">
      <div className="value">{typeof value === "number" ? value.toLocaleString() : value}</div>
      <div className="label">{label}</div>
    </div>
  );
  return (
    <div className="overview">
      {stat("variants", o.variants)}
      {stat("ClinVar", o.clinvar_annotated)}
      {stat("gnomAD", o.gnomad_annotated)}
      {stat("PGx genes", o.pgx_genes)}
      {stat("CNV", o.cnv)}
      {stat("SV", o.sv)}
      {stat("build", o.meta.source_vcf ? (o.meta.build ?? "GRCh38") : "—")}
    </div>
  );
}

function useAsync<T>() {
  const [data, setData] = useState<T | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const run = async (fn: () => Promise<T>) => {
    setLoading(true);
    setError(null);
    try {
      setData(await fn());
    } catch (e) {
      setError(String((e as Error).message ?? e));
      setData(null);
    } finally {
      setLoading(false);
    }
  };
  return { data, loading, error, run };
}

function SearchView() {
  const [q, setQ] = useState("");
  const { data, loading, error, run } = useAsync<VariantPage>();

  const search = () => {
    const kind = classifyQuery(q);
    run(() => (kind === "rsid" ? api.byRsid(q) : kind === "region" ? api.byRegion(q) : api.byGene(q)));
  };

  return (
    <section>
      <div className="searchbar">
        <input
          value={q}
          placeholder="gene (BRCA1), rsID (rs1799853), or region (chr7:117480000-117670000)"
          onChange={(e) => setQ(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && search()}
        />
        <button onClick={search}>Search</button>
      </div>
      <p className="hint">Detected: <code>{q ? classifyQuery(q) : "—"}</code></p>
      {loading && <p>Searching…</p>}
      {error && <div className="banner error">{error}</div>}
      {data && <VariantTable page={data} />}
    </section>
  );
}

function ClinicalView() {
  const [gene, setGene] = useState("");
  const [sig, setSig] = useState("");
  const { data, loading, error, run } = useAsync<VariantPage>();
  useEffect(() => {
    run(() => api.clinical());
  }, []);
  return (
    <section>
      <div className="searchbar">
        <input value={gene} placeholder="gene (optional)" onChange={(e) => setGene(e.target.value)} />
        <input value={sig} placeholder="significance (default: pathogenic)" onChange={(e) => setSig(e.target.value)} />
        <button onClick={() => run(() => api.clinical(gene, sig))}>Filter</button>
      </div>
      {loading && <p>Loading…</p>}
      {error && <div className="banner error">{error}</div>}
      {data && <VariantTable page={data} clinical />}
    </section>
  );
}

function PgxView() {
  const [gene, setGene] = useState("");
  const [drug, setDrug] = useState("");
  const { data, loading, error, run } = useAsync<PgxResult>();
  useEffect(() => {
    run(() => api.pgx());
  }, []);
  return (
    <section>
      <div className="searchbar">
        <input value={gene} placeholder="gene (CYP2C19)" onChange={(e) => setGene(e.target.value)} />
        <input value={drug} placeholder="drug (clopidogrel)" onChange={(e) => setDrug(e.target.value)} />
        <button onClick={() => run(() => api.pgx(gene, drug))}>Filter</button>
      </div>
      {loading && <p>Loading…</p>}
      {error && <div className="banner error">{error}</div>}
      {data && (
        <>
          <h3>Gene diplotypes</h3>
          {data.genes.length === 0 ? (
            <p className="hint">No PGx results. Run the PharmCAT annotation step.</p>
          ) : (
            <table>
              <thead>
                <tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th><th>Activity</th></tr>
              </thead>
              <tbody>
                {data.genes.map((g, i) => (
                  <tr key={i}><td>{g.gene}</td><td>{g.diplotype}</td><td>{g.phenotype}</td><td>{g.activity_score}</td></tr>
                ))}
              </tbody>
            </table>
          )}
          {data.drugs.length > 0 && (
            <>
              <h3>Drug guidance</h3>
              <table>
                <thead>
                  <tr><th>Drug</th><th>Gene</th><th>Source</th><th>Recommendation</th></tr>
                </thead>
                <tbody>
                  {data.drugs.map((d, i) => (
                    <tr key={i}><td>{d.drug}</td><td>{d.gene}</td><td>{d.source}</td><td>{d.recommendation}</td></tr>
                  ))}
                </tbody>
              </table>
            </>
          )}
        </>
      )}
    </section>
  );
}

function SqlView() {
  const [q, setQ] = useState("SELECT chrom, pos, ref, alt, rsid, gene, clnsig FROM variants LIMIT 50");
  const { data, loading, error, run } = useAsync<SqlResult>();
  return (
    <section>
      <textarea value={q} onChange={(e) => setQ(e.target.value)} rows={4} />
      <button onClick={() => run(() => api.sql(q))}>Run (read-only)</button>
      {loading && <p>Running…</p>}
      {error && <div className="banner error">{error}</div>}
      {data && (
        <table>
          <thead>
            <tr>{data.columns.map((c) => <th key={c}>{c}</th>)}</tr>
          </thead>
          <tbody>
            {data.rows.map((row, i) => (
              <tr key={i}>{row.map((cell, j) => <td key={j}>{cell === null ? "" : String(cell)}</td>)}</tr>
            ))}
          </tbody>
        </table>
      )}
      {data && <p className="hint">{data.rows.length} rows (capped at {data.truncated_to}).</p>}
    </section>
  );
}

function VariantTable({ page, clinical }: { page: VariantPage; clinical?: boolean }) {
  if (page.hits.length === 0) return <p className="hint">No matches ({page.total} total).</p>;
  return (
    <>
      <p className="hint">
        Showing {page.hits.length} of {page.total.toLocaleString()} matches.
      </p>
      <table>
        <thead>
          <tr>
            <th>Locus</th><th>Ref→Alt</th><th>rsID</th><th>GT</th><th>Gene</th>
            {clinical ? <th>Significance</th> : <th>Consequence</th>}
            {clinical ? <th>Disease</th> : <th>gnomAD AF</th>}
          </tr>
        </thead>
        <tbody>
          {page.hits.map((v: Variant, i) => (
            <tr key={i}>
              <td className="mono">{v.chrom}:{v.pos.toLocaleString()}</td>
              <td className="mono">{v.ref}→{v.alt}</td>
              <td>{v.rsid ?? ""}</td>
              <td className="mono">{v.gt ?? ""}</td>
              <td>{v.gene ?? ""}</td>
              {clinical ? (
                <td className={/pathogenic/i.test(v.clnsig ?? "") ? "sig-path" : ""}>{v.clnsig ?? ""}</td>
              ) : (
                <td>{v.consequence ?? ""}</td>
              )}
              {clinical ? <td>{v.clndn ?? ""}</td> : <td>{v.gnomad_af ?? ""}</td>}
            </tr>
          ))}
        </tbody>
      </table>
    </>
  );
}
