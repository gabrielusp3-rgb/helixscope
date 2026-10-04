"use client";

import { useMemo, useState } from "react";
import { formatScientificValue } from "@/lib/api/numeric";
import { downloadCsv, downloadText } from "@/lib/export/download";

const WINDOW = 250;

export function DataTable({
  rows,
  caption,
  testId,
  emptyLabel = "0 rows.",
  onRowSelect,
  selectedKey,
}: {
  rows: Array<Record<string, unknown>>;
  caption?: string;
  testId?: string;
  emptyLabel?: string;
  onRowSelect?: (row: Record<string, unknown>, index: number) => void;
  selectedKey?: string | null;
}) {
  const columns = useMemo(() => {
    const keys = new Set<string>();
    rows.forEach((row) => Object.keys(row).forEach((key) => keys.add(key)));
    return Array.from(keys);
  }, [rows]);
  const [sortKey, setSortKey] = useState<string | null>(null);
  const [dir, setDir] = useState<1 | -1>(1);
  const [query, setQuery] = useState("");
  const sorted = useMemo(() => {
    const needle = query.trim().toLowerCase();
    const filtered = needle
      ? rows.filter((row) =>
          columns.some((col) => formatScientificValue(row[col], "").toLowerCase().includes(needle)),
        )
      : rows;
    if (!sortKey) return filtered;
    return [...filtered].sort((a, b) => {
      const av = String(a[sortKey] ?? "");
      const bv = String(b[sortKey] ?? "");
      return av.localeCompare(bv, undefined, { numeric: true }) * dir;
    });
  }, [rows, sortKey, dir, query, columns]);
  const truncated = sorted.length > WINDOW;
  const visible = truncated ? sorted.slice(0, WINDOW) : sorted;

  if (rows.length === 0) {
    return (
      <p className="empty-state" data-testid={testId}>
        {emptyLabel}
      </p>
    );
  }

  return (
    <div className="hs-table-wrap" data-testid={testId}>
      {caption ? <p className="hs-topbar-meta">{caption}</p> : null}
      <div className="hs-table-toolbar">
        <label className="hs-table-search">
          <span className="sr-only">Filter table</span>
          <input
            type="search"
            value={query}
            onChange={(event) => setQuery(event.target.value)}
            placeholder="Filter rows"
            aria-label={caption ? `Filter ${caption}` : "Filter table"}
          />
        </label>
        <button
          type="button"
          className="btn-ghost"
          onClick={() => {
            const header = columns.join("\t");
            const body = sorted
              .map((row) => columns.map((col) => formatScientificValue(row[col], "")).join("\t"))
              .join("\n");
            void navigator.clipboard?.writeText(`${header}\n${body}\n`);
          }}
        >
          Copy TSV
        </button>
        <button
          type="button"
          className="btn-ghost"
          onClick={() => downloadCsv(`${(caption || "table").replace(/\s+/g, "_")}.csv`, sorted)}
        >
          Export CSV
        </button>
        <button
          type="button"
          className="btn-ghost"
          onClick={() =>
            downloadText(
              `${(caption || "table").replace(/\s+/g, "_")}.tsv`,
              `${columns.join("\t")}\n${sorted
                .map((row) => columns.map((col) => formatScientificValue(row[col], "")).join("\t"))
                .join("\n")}\n`,
              "text/tab-separated-values;charset=utf-8",
            )
          }
        >
          Export TSV
        </button>
      </div>
      {truncated ? (
        <p className="hs-topbar-meta">
          Showing {visible.length} of {sorted.length} rows. Sorting and filtering are presentational.
          Download CSV/TSV for the full API table. Values are not recomputed.
        </p>
      ) : null}
      <div className="hs-table-scroll">
        <table className="hs-table">
          {caption ? <caption className="sr-only">{caption}</caption> : null}
          <thead>
            <tr>
              {columns.map((col) => (
                <th key={col} scope="col">
                  <button
                    type="button"
                    className="btn-ghost"
                    onClick={() => {
                      if (sortKey === col) setDir((d) => (d === 1 ? -1 : 1));
                      else {
                        setSortKey(col);
                        setDir(1);
                      }
                    }}
                  >
                    {col}
                  </button>
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {visible.map((row, index) => (
              <tr
                key={index}
                data-clickable={onRowSelect ? "true" : "false"}
                data-selected={selectedKey != null && String(row.start ?? row.identifier ?? index) === selectedKey ? "true" : "false"}
                onClick={() => onRowSelect?.(row, index)}
              >
                {columns.map((col) => (
                  <td key={col}>{formatScientificValue(row[col], "")}</td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
