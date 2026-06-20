import { geoMercator, geoPath } from "d3-geo";
import { useEffect, useState } from "react";
import { feature } from "topojson-client";
import worldData from "world-atlas/countries-110m.json";
import {
  api,
  type AncestryComponent,
  type AncestrySummary,
  type AssociationPage,
  classifyQuery,
  type Overview,
  type PcaPoint,
  type PgsResult,
  type PgxResult,
  type SqlResult,
  type TraitsReport,
  type Variant,
  type VariantPage,
  type WhatsNew,
} from "./api";

type Tab = "search" | "clinical" | "pgx" | "ancestry" | "risk" | "traits" | "gwas" | "changelog" | "sql";

// Continental colors for the ancestry PCA cloud.
const SUPERPOP_COLORS: Record<string, string> = {
  EUR: "#4ea1ff", AFR: "#ffcc4e", EAS: "#7c6cff", SAS: "#4ecb8b", AMR: "#ff8a4e",
  MID: "#e36bd0", OCE: "#2fd0c7",
};
const SUPERPOP_NAMES: Record<string, string> = {
  EUR: "European", AFR: "African", EAS: "East Asian", SAS: "South Asian", AMR: "American",
  MID: "Middle Eastern", OCE: "Oceanian",
};

// Approximate ancestral-homeland coordinates [lng, lat] + continent for reference populations,
// so the geographic map can place a marker per population the genome is closest to. (1000 Genomes
// codes use the ancestral homeland, not the sampling site — e.g. CEU is placed in NW Europe.)
const POP_COORDS: Record<string, { lng: number; lat: number; group: string }> = {
  // 1000 Genomes — European
  GBR: { lng: -1.5, lat: 53, group: "EUR" }, CEU: { lng: 7, lat: 51, group: "EUR" },
  FIN: { lng: 25, lat: 62, group: "EUR" }, IBS: { lng: -3.7, lat: 40, group: "EUR" },
  TSI: { lng: 11, lat: 43, group: "EUR" },
  // HGDP — European
  French: { lng: 2.3, lat: 46.5, group: "EUR" }, Orcadian: { lng: -3, lat: 59, group: "EUR" },
  Basque: { lng: -2, lat: 43, group: "EUR" }, Sardinian: { lng: 9, lat: 40, group: "EUR" },
  BergamoItalian: { lng: 9.7, lat: 45.7, group: "EUR" }, Tuscan: { lng: 11, lat: 43.3, group: "EUR" },
  Russian: { lng: 37, lat: 56, group: "EUR" }, Adygei: { lng: 40, lat: 44, group: "MID" },
  // 1000 Genomes — African / East Asian / South Asian / American
  YRI: { lng: 8, lat: 7, group: "AFR" }, LWK: { lng: 37, lat: 0, group: "AFR" },
  GWD: { lng: -15, lat: 13, group: "AFR" }, MSL: { lng: -12, lat: 8, group: "AFR" },
  ESN: { lng: 8, lat: 9, group: "AFR" }, ACB: { lng: -59, lat: 13, group: "AFR" },
  ASW: { lng: -98, lat: 38, group: "AFR" },
  CHB: { lng: 116, lat: 40, group: "EAS" }, JPT: { lng: 139, lat: 36, group: "EAS" },
  CHS: { lng: 113, lat: 28, group: "EAS" }, CDX: { lng: 100, lat: 22, group: "EAS" },
  KHV: { lng: 106, lat: 11, group: "EAS" },
  GIH: { lng: 72, lat: 23, group: "SAS" }, PJL: { lng: 74, lat: 31, group: "SAS" },
  BEB: { lng: 90, lat: 24, group: "SAS" }, STU: { lng: 81, lat: 8, group: "SAS" },
  ITU: { lng: 79, lat: 13, group: "SAS" },
  MXL: { lng: -102, lat: 23, group: "AMR" }, PUR: { lng: -66, lat: 18, group: "AMR" },
  CLM: { lng: -74, lat: 4, group: "AMR" }, PEL: { lng: -77, lat: -12, group: "AMR" },
};

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
        {(["search", "clinical", "pgx", "ancestry", "risk", "traits", "gwas", "changelog", "sql"] as Tab[]).map((t) => (
          <button key={t} className={tab === t ? "active" : ""} onClick={() => setTab(t)}>
            {t === "pgx" ? "Pharmacogenomics" : t[0].toUpperCase() + t.slice(1)}
          </button>
        ))}
      </nav>

      <main>
        {tab === "search" && <SearchView />}
        {tab === "clinical" && <ClinicalView />}
        {tab === "pgx" && <PgxView />}
        {tab === "ancestry" && <AncestryView />}
        {tab === "risk" && <RiskView />}
        {tab === "traits" && <TraitsView />}
        {tab === "gwas" && <GwasView />}
        {tab === "changelog" && <ChangelogView />}
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

function AncestryMap({ populations }: { populations: AncestryComponent[] }) {
  const W = 720, H = 420;
  /* eslint-disable @typescript-eslint/no-explicit-any */
  const world: any = feature(worldData as any, (worldData as any).objects.countries);
  const cand = populations
    .map((p) => ({ p, c: POP_COORDS[p.code] }))
    .filter((m): m is { p: AncestryComponent; c: { lng: number; lat: number; group: string } } => !!m.c);
  // Zoom to the markers' region: fit the projection to a MultiPoint of the homelands, enforcing a
  // minimum span so a tight cluster still shows surrounding context. Falls back to a broad world
  // view if no homelands are known.
  const coords: [number, number][] = cand.map((m) => [m.c.lng, m.c.lat]);
  if (!coords.length) {
    coords.push([-150, -45], [160, 70]);
  } else {
    const lons = coords.map((c) => c[0]), lats = coords.map((c) => c[1]);
    const cLon = (Math.min(...lons) + Math.max(...lons)) / 2;
    const cLat = (Math.min(...lats) + Math.max(...lats)) / 2;
    const MIN = 18;  // degrees
    if (Math.max(Math.max(...lons) - Math.min(...lons), Math.max(...lats) - Math.min(...lats)) < MIN) {
      coords.push([cLon - MIN / 2, cLat - MIN / 2], [cLon + MIN / 2, cLat + MIN / 2]);
    }
  }
  const proj = geoMercator()
    .fitExtent([[150, 90], [W - 150, H - 90]], { type: "MultiPoint", coordinates: coords } as any)
    .clipExtent([[0, 0], [W, H]]);  // clip far geometry so Mercator's poles can't blow up
  const path = geoPath(proj);
  /* eslint-enable @typescript-eslint/no-explicit-any */
  const markers = cand.map((m) => ({ ...m, xy: proj([m.c.lng, m.c.lat]) })).filter((m) => m.xy);
  return (
    <>
      <h3 style={{ marginTop: "1.5rem" }}>Ancestral homelands (geographic)</h3>
      <svg viewBox={`0 0 ${W} ${H}`} preserveAspectRatio="xMidYMid meet" className="geo">
        {world.features.map((f: unknown, i: number) => (
          <path key={i} d={path(f as never) ?? ""} className="geo-land" />
        ))}
        {markers.map((m) => {
          const [x, y] = m.xy as [number, number];
          const r = 6 + m.p.proportion * 38;
          return (
            <g key={m.p.code}>
              <circle cx={x} cy={y} r={r} fill={SUPERPOP_COLORS[m.c.group] ?? "#fff"}
                opacity={0.5} stroke="white" strokeWidth={1.2} />
              <text x={x + r + 4} y={y + 3} className="geo-lbl">
                {m.p.name.split(" (")[0]} {(m.p.proportion * 100).toFixed(0)}%
              </text>
            </g>
          );
        })}
      </svg>
      <p className="hint">Auto-zoomed to where your ancestry sits. Markers are at each population's
        ancestral homeland, sized by your sub-continental proportion — approximate, reflecting genetic
        similarity, not a precise birthplace.</p>
    </>
  );
}

function AncestryView() {
  const { data, loading, error, run } = useAsync<AncestrySummary>();
  useEffect(() => {
    run(() => api.ancestry());
  }, []);
  if (loading) return <p>Loading…</p>;
  if (error) return <div className="banner error">{error}</div>;
  if (!data || data.components.length === 0)
    return <p className="hint">No ancestry results yet. Run <code>locus ancestry</code>.</p>;

  // PCA scatter scaling
  const pts = data.pca;
  const xs = pts.map((p) => p.pc1);
  const ys = pts.map((p) => p.pc2);
  const [xmin, xmax] = [Math.min(...xs), Math.max(...xs)];
  const [ymin, ymax] = [Math.min(...ys), Math.max(...ys)];
  const W = 460, H = 320, pad = 34;
  const sx = (v: number) => pad + ((v - xmin) / (xmax - xmin || 1)) * (W - 2 * pad);
  const sy = (v: number) => H - pad - ((v - ymin) / (ymax - ymin || 1)) * (H - 2 * pad);

  // Continent halos + labels: group reference points by continent and draw a soft blob behind each.
  const groups: Record<string, PcaPoint[]> = {};
  for (const p of pts.filter((q) => !q.is_sample)) (groups[p.group ?? "?"] ??= []).push(p);
  const halos = Object.entries(groups).map(([g, ps]) => {
    const cx = ps.reduce((s, p) => s + sx(p.pc1), 0) / ps.length;
    const cy = ps.reduce((s, p) => s + sy(p.pc2), 0) / ps.length;
    const rms = Math.sqrt(ps.reduce((s, p) => s + (sx(p.pc1) - cx) ** 2 + (sy(p.pc2) - cy) ** 2, 0) / ps.length);
    return { g, cx, cy, r: Math.max(20, rms * 1.6 + 12), color: SUPERPOP_COLORS[g] ?? "#8b92a7" };
  });
  const x0 = sx(0), y0 = sy(0);  // origin axes (PC space is centered near 0)

  return (
    <section>
      <h3>Continental ancestry</h3>
      {data.components.map((c) => (
        <div key={c.code} className="bar-row">
          <span className="bar-label">{c.name}</span>
          <span className="bar-track"><span className="bar-fill" style={{ width: `${c.proportion * 100}%` }} /></span>
          <span className="bar-val">{(c.proportion * 100).toFixed(0)}%</span>
        </div>
      ))}
      <h3>Sub-continental (closest populations)</h3>
      {data.populations.map((c) => (
        <div key={c.code} className="bar-row">
          <span className="bar-label">{c.name}</span>
          <span className="bar-track"><span className="bar-fill" style={{ width: `${c.proportion * 100}%` }} /></span>
          <span className="bar-val">{(c.proportion * 100).toFixed(0)}%</span>
        </div>
      ))}
      <AncestryMap populations={data.populations} />

      <h3 style={{ marginTop: "1.5rem" }}>Genetic-similarity map (PC1 × PC2)</h3>
      <svg width={W} height={H} className="pca">
        {/* soft continent regions */}
        {halos.map((h) => (
          <circle key={`halo-${h.g}`} cx={h.cx} cy={h.cy} r={h.r} fill={h.color} opacity={0.1} />
        ))}
        {/* origin axes */}
        <line x1={x0} y1={pad} x2={x0} y2={H - pad} className="pca-axis" />
        <line x1={pad} y1={y0} x2={W - pad} y2={y0} className="pca-axis" />
        <text x={W - pad} y={y0 - 5} className="pca-axislbl" textAnchor="end">PC1 →</text>
        <text x={x0 + 5} y={pad + 2} className="pca-axislbl">↑ PC2</text>
        {/* continent labels */}
        {halos.map((h) => (
          <text key={`lbl-${h.g}`} x={h.cx} y={h.cy - h.r - 2} className="pca-grouplbl"
            fill={h.color} textAnchor="middle">{h.g}</text>
        ))}
        {pts.filter((p) => !p.is_sample).map((p) => (
          <circle key={p.label} cx={sx(p.pc1)} cy={sy(p.pc2)} r={5}
            fill={SUPERPOP_COLORS[p.group ?? ""] ?? "#8b92a7"} opacity={0.85}>
            <title>{p.label} ({p.group})</title>
          </circle>
        ))}
        {pts.filter((p) => p.is_sample).map((p) => (
          <g key="you">
            <circle cx={sx(p.pc1)} cy={sy(p.pc2)} r={7} className="pca-you" />
            <text x={sx(p.pc1) + 9} y={sy(p.pc2) + 4} className="pca-lbl you">you</text>
          </g>
        ))}
      </svg>
      <div className="pca-legend">
        {Object.entries(SUPERPOP_COLORS).map(([code, color]) => (
          <span key={code} className="legend-item">
            <span className="legend-dot" style={{ background: color }} /> {SUPERPOP_NAMES[code] ?? code}
          </span>
        ))}
        <span className="legend-item"><span className="legend-dot" style={{ background: "var(--path)" }} /> you</span>
      </div>
      <p className="hint">Each dot is a reference population (1000 Genomes + HGDP), colored by continent;
        your genome (red) sits among them. This is genetic-similarity space (PC1×PC2), not geography.
        {" "}{data.note}</p>
    </section>
  );
}

function RiskView() {
  const { data, loading, error, run } = useAsync<PgsResult[]>();
  useEffect(() => {
    run(() => api.pgs());
  }, []);
  if (loading) return <p>Loading…</p>;
  if (error) return <div className="banner error">{error}</div>;
  if (!data || data.length === 0)
    return <p className="hint">No polygenic scores yet. Run <code>locus ancestry</code>.</p>;
  return (
    <section>
      <h3>Polygenic (aggregate) risk</h3>
      {data.map((s) => (
        <div key={s.pgs_id} className="bar-row">
          <span className="bar-label">{s.trait}</span>
          <span className="bar-track">
            {s.percentile !== null && <span className="bar-fill" style={{ width: `${s.percentile}%` }} />}
          </span>
          <span className="bar-val">
            {s.percentile !== null ? `${s.percentile.toFixed(0)}th pct` : "raw only"}
          </span>
        </div>
      ))}
      <p className="hint">
        Percentiles are within your ancestry-matched 1000 Genomes reference{data[0]?.ancestry ? ` (${data[0].ancestry})` : ""};
        they're research-grade estimates, not diagnoses, and absolute risk across ancestries is unreliable.
        Coverage shows the fraction of each score's variants callable in your genome.
      </p>
    </section>
  );
}

function GwasView() {
  const [q, setQ] = useState("type 2 diabetes");
  const { data, loading, error, run } = useAsync<AssociationPage>();
  useEffect(() => {
    run(() => api.gwas(q));
  }, []);
  return (
    <section>
      <div className="search-row">
        <input value={q} onChange={(e) => setQ(e.target.value)} placeholder="trait, e.g. height, asthma…"
          onKeyDown={(e) => e.key === "Enter" && run(() => api.gwas(q))} />
        <button onClick={() => run(() => api.gwas(q))}>Search</button>
      </div>
      {loading && <p>Loading…</p>}
      {error && <div className="banner error">{error}</div>}
      {data && (
        <>
          <p className="hint">{data.total.toLocaleString()} carried risk allele(s) for "{data.trait}". {data.note}</p>
          <table>
            <thead>
              <tr><th>rsID</th><th>Risk</th><th>Zygosity</th><th>Trait</th><th>OR/β</th><th>P</th></tr>
            </thead>
            <tbody>
              {data.hits.map((a, i) => (
                <tr key={i}>
                  <td>{a.rsid}</td>
                  <td className="mono">{a.risk_allele}</td>
                  <td>{a.zygosity}</td>
                  <td>{a.mapped_trait}</td>
                  <td className="mono">{a.or_beta ?? ""}</td>
                  <td className="mono">{a.pval.toExponential(0)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </>
      )}
    </section>
  );
}

function TraitsView() {
  const traits = useAsync<TraitsReport>();
  const acmg = useAsync<VariantPage>();
  useEffect(() => {
    traits.run(() => api.traits());
    acmg.run(() => api.secondaryFindings());
  }, []);
  return (
    <section>
      <h3>ACMG secondary findings</h3>
      {acmg.loading && <p>Loading…</p>}
      {acmg.error && <div className="banner error">{acmg.error}</div>}
      {acmg.data && acmg.data.total === 0 && (
        <p className="hint">None — no pathogenic/likely-pathogenic variants in the ACMG SF v3.2
          actionable genes. A reassuring result (confirm clinically if ever flagged).</p>
      )}
      {acmg.data && acmg.data.total > 0 && <VariantTable page={acmg.data} clinical />}

      <h3 style={{ marginTop: "1.5rem" }}>Traits & wellness</h3>
      {traits.loading && <p>Loading…</p>}
      {traits.error && <div className="banner error">{traits.error}</div>}
      {traits.data && traits.data.total === 0 && (
        <p className="hint">No traits computed yet. Run <code>locus traits</code>.</p>
      )}
      {traits.data && traits.data.total > 0 && (
        <table>
          <thead>
            <tr><th>Trait</th><th>Genotype</th><th>Interpretation</th><th>rsID</th></tr>
          </thead>
          <tbody>
            {traits.data.traits.map((t) => (
              <tr key={t.rsid}>
                <td>{t.trait}{t.category === "pharmacogenomic" ? " ⚕" : ""}</td>
                <td className="mono">{t.genotype}</td>
                <td>{t.interpretation}</td>
                <td>{t.rsid}</td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
      {traits.data && traits.data.total > 0 && <p className="hint">{traits.data.note}</p>}
    </section>
  );
}

function ChangelogView() {
  const { data, loading, error, run } = useAsync<WhatsNew>();
  useEffect(() => {
    run(() => api.whatsNew());
  }, []);
  if (loading) return <p>Loading…</p>;
  if (error) return <div className="banner error">{error}</div>;
  if (!data || data.total === 0)
    return <p className="hint">Nothing new yet. Run <code>locus refresh</code> to check for newly-published findings (ClinVar reclassifications, new polygenic scores).</p>;
  const order = ["strong", "moderate", "weak", "info"];
  return (
    <section>
      <h3>What's new in your genome</h3>
      <p className="hint">
        {data.total} finding(s) from the last <code>locus refresh</code>
        {Object.entries(data.counts_by_tier).map(([t, n]) => ` · ${n} ${t}`).join("")}
      </p>
      {order.filter((t) => data.findings.some((f) => f.tier === t)).map((tier) => (
        <div key={tier}>
          <h4 className={tier === "strong" ? "sig-path" : ""}>{tier[0].toUpperCase() + tier.slice(1)}</h4>
          {data.findings.filter((f) => f.tier === tier).map((f, i) => (
            <div key={i} className="bar-row">
              <span className="bar-label">{f.title}</span>
              <span className="bar-val">
                {f.chrom ? `${f.chrom}:${f.pos?.toLocaleString()}` : f.source}
              </span>
              {f.detail && <span className="hint" style={{ flexBasis: "100%" }}>{f.detail}</span>}
            </div>
          ))}
        </div>
      ))}
      <p className="hint">{data.note}</p>
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
            <th>AlphaMissense</th>
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
              <td className={/pathogenic/.test(v.am_class ?? "") ? "sig-path" : ""}>
                {v.am_class ? `${v.am_class}${v.am_pathogenicity != null ? ` (${v.am_pathogenicity.toFixed(2)})` : ""}` : ""}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </>
  );
}
