// Typed client for the Locus FastAPI backend.

export interface Variant {
  chrom: string;
  pos: number;
  ref: string;
  alt: string;
  rsid: string | null;
  gt: string | null;
  filter: string | null;
  gene: string | null;
  consequence: string | null;
  clnsig: string | null;
  clndn: string | null;
  clnrevstat: string | null;
  gnomad_af: number | null;
  gnomad_af_grpmax: number | null;
}

export interface VariantPage {
  total: number;
  limit: number;
  offset: number;
  hits: Variant[];
}

export interface PgxGene {
  gene: string;
  diplotype: string | null;
  phenotype: string | null;
  activity_score: string | null;
}
export interface PgxDrug {
  drug: string;
  gene: string | null;
  source: string | null;
  recommendation: string | null;
}
export interface PgxResult {
  genes: PgxGene[];
  drugs: PgxDrug[];
}

export interface Overview {
  meta: Record<string, string>;
  variants: number;
  clinvar_annotated: number;
  gnomad_annotated: number;
  pgx_genes: number;
  cnv: number;
  sv: number;
}

export interface SqlResult {
  columns: string[];
  rows: unknown[][];
  truncated_to: number;
}

async function get<T>(path: string): Promise<T> {
  const r = await fetch(path);
  if (!r.ok) throw new Error((await r.json().catch(() => ({ detail: r.statusText }))).detail ?? r.statusText);
  return r.json();
}

export const api = {
  overview: () => get<Overview>("/api/overview"),
  byRsid: (rsid: string) => get<VariantPage>(`/api/variant/rsid/${encodeURIComponent(rsid)}`),
  byGene: (gene: string) => get<VariantPage>(`/api/gene/${encodeURIComponent(gene)}`),
  byRegion: (region: string) => get<VariantPage>(`/api/region?region=${encodeURIComponent(region)}`),
  clinical: (gene = "", significance = "") =>
    get<VariantPage>(`/api/clinical?gene=${encodeURIComponent(gene)}&significance=${encodeURIComponent(significance)}`),
  pgx: (gene = "", drug = "") =>
    get<PgxResult>(`/api/pgx?gene=${encodeURIComponent(gene)}&drug=${encodeURIComponent(drug)}`),
  sql: async (query: string): Promise<SqlResult> => {
    const r = await fetch("/api/sql", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ query }),
    });
    if (!r.ok) throw new Error((await r.json().catch(() => ({ detail: r.statusText }))).detail ?? r.statusText);
    return r.json();
  },
};

// Heuristic: decide which lookup to run from a free-text query.
export function classifyQuery(q: string): "rsid" | "region" | "gene" {
  const s = q.trim();
  if (/^rs\d+$/i.test(s)) return "rsid";
  if (/^(chr)?[\w]+:[\d,]+(-[\d,]+)?$/i.test(s)) return "region";
  return "gene";
}
