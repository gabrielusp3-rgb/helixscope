"use client";

import { useMemo, useState } from "react";
import { asList, asRecord } from "@/lib/api/numeric";
import { residueColor } from "@/lib/visual/nucleotides";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";

const COL_WINDOW = 240;
const ROW_WINDOW = 80;

export function MsaViewer({ result }: { result: Record<string, unknown> }) {
  const rows = asList(result.rows).map((row) => asRecord(row));
  const conservation = asRecord(result.conservation);
  const scores = asList(conservation.scores).map((v) => Number(v));
  const consensus = asRecord(result.consensus);
  const { selectedMsaId, setSelectedMsaId } = useWorkspace();
  const consSeq = String(consensus.sequence || "");
  const length = String(asRecord(rows[0]).aligned || consSeq).length;
  const molecule = String(result.molecule || "DNA");
  const selectedCol = Number.parseInt(selectedMsaId ?? "", 10);
  const [font, setFont] = useState(12);
  const [query, setQuery] = useState("");
  const selected = Number.isFinite(selectedCol) ? selectedCol : 0;
  const needsWindow = length > COL_WINDOW;
  const start = needsWindow ? Math.max(0, Math.min(length - COL_WINDOW, selected - Math.floor(COL_WINDOW / 2))) : 0;
  const end = needsWindow ? Math.min(length, start + COL_WINDOW) : length;
  const filtered = useMemo(() => {
    const needle = query.trim().toLowerCase();
    if (!needle) return rows;
    return rows.filter((row) => String(row.identifier || "").toLowerCase().includes(needle));
  }, [rows, query]);
  const truncatedRows = filtered.length > ROW_WINDOW;
  const visibleRows = truncatedRows ? filtered.slice(0, ROW_WINDOW) : filtered;

  function sliceSeq(seq: string): string {
    return needsWindow ? seq.slice(start, end) : seq;
  }

  function residueSpan(seq: string, keyPrefix: string) {
    return sliceSeq(seq)
      .split("")
      .map((ch, offset) => {
        const i = start + offset;
        const gap = ch === "-" || ch === ".";
        return (
          <span
            key={`${keyPrefix}-${i}`}
            data-selected={i === selectedCol ? "true" : "false"}
            data-gap={gap ? "true" : "false"}
            style={{ color: residueColor(ch, molecule), fontSize: font }}
            onClick={() => setSelectedMsaId(String(i))}
          >
            {ch}
          </span>
        );
      });
  }

  const ruler = sliceSeq(" ".repeat(length) || consSeq)
    .split("")
    .map((_, offset) => {
      const i = start + offset + 1;
      return (
        <span key={i} style={{ fontSize: Math.max(8, font - 3), color: i % 10 === 0 ? "#AFCBF1" : "#3d4d66" }}>
          {i % 10 === 0 ? String(i % 100).padStart(2, "0").slice(-1) : "·"}
        </span>
      );
    });

  return (
    <div className="msa-viewer" data-testid="msa-viewer">
      <div className="msa-toolbar">
        <label>
          Font
          <input
            type="range"
            min={8}
            max={18}
            value={font}
            onChange={(event) => setFont(Number(event.target.value))}
            aria-label="MSA font size"
          />
        </label>
        <label>
          Search sequence ID
          <input
            type="search"
            value={query}
            onChange={(event) => setQuery(event.target.value)}
            data-testid="msa-id-search"
            placeholder="identifier"
          />
        </label>
        <p className="hs-topbar-meta">
          {needsWindow ? `columns ${start + 1}-${end} of ${length}` : `${length} columns`} · {filtered.length} sequences
        </p>
      </div>
      <div className="msa-row">
        <div className="msa-id">ruler</div>
        <div className="msa-ruler">{ruler}</div>
      </div>
      <div className="msa-row">
        <div className="msa-id">consensus</div>
        <div className="msa-seq">{residueSpan(consSeq, "cons")}</div>
      </div>
      {visibleRows.map((row) => (
        <div className="msa-row" key={String(row.identifier)}>
          <div className="msa-id">{String(row.identifier)}</div>
          <div className="msa-seq">{residueSpan(String(row.aligned || ""), String(row.identifier))}</div>
        </div>
      ))}
      {truncatedRows ? (
        <p className="hs-topbar-meta" style={{ padding: 8 }}>
          Showing {visibleRows.length} of {filtered.length} sequences (virtualized). Filter by ID to inspect others.
        </p>
      ) : null}
      <div className="hs-topbar-meta" style={{ padding: 8 }}>
        Conservation scores come from the API ({scores.length} values). Selecting a column is local UI
        state and does not recompute conservation. Residue colors are presentation only.
      </div>
    </div>
  );
}
