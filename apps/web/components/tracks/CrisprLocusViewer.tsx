"use client";

import { asRecord, chartNumeric, formatScientificValue } from "@/lib/api/numeric";
import { residueColor } from "@/lib/visual/nucleotides";

export function CrisprLocusViewer({
  sequence,
  guides,
  selectedIndex,
  onSelect,
  pamLabel,
  cutLabel,
}: {
  sequence: string;
  guides: Array<Record<string, unknown>>;
  selectedIndex: number;
  onSelect: (index: number) => void;
  pamLabel: string;
  cutLabel: string;
}) {
  const length = sequence.length || 1;
  const width = 920;
  const height = 120;
  const x = (pos: number) => 16 + (pos / length) * (width - 32);
  const selected = guides[selectedIndex] ? asRecord(guides[selectedIndex]) : {};
  const cut = chartNumeric(selected.cut_site ?? selected.nick_site ?? selected.edit_window_start);

  return (
    <section className="hs-panel science-surface" data-testid="crispr-locus">
      <h2>Locus</h2>
      <p className="hs-topbar-meta">
        {pamLabel}. {cutLabel}. Positions are from the API result, not recomputed in the browser.
      </p>
      <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label="CRISPR locus from API coordinates">
        <rect width={width} height={height} fill="#01040B" />
        <line x1="16" x2={width - 16} y1="36" y2="36" stroke="rgba(175,203,241,0.4)" />
        {guides.map((row, index) => {
          const guide = asRecord(row);
          const start = chartNumeric(guide.start_0based ?? guide.start ?? guide.spacer_start);
          const end = chartNumeric(guide.end_0based ?? guide.end ?? guide.spacer_end);
          if (start === null || end === null) return null;
          const active = index === selectedIndex;
          return (
            <rect
              key={index}
              x={x(Math.min(start, end))}
              y={active ? 48 : 54}
              width={Math.max(3, Math.abs(x(end) - x(start)))}
              height={active ? 18 : 10}
              rx="3"
              fill={active ? "#D9E4FC" : "#2355A0"}
              tabIndex={0}
              role="button"
              aria-label={`Guide ${index + 1}`}
              onClick={() => onSelect(index)}
            />
          );
        })}
        {cut !== null ? <line x1={x(cut)} x2={x(cut)} y1="20" y2="90" stroke="#ff6b8a" strokeDasharray="3 3" /> : null}
        {length <= 80
          ? sequence.split("").map((ch, index) => (
              <text key={index} x={x(index)} y="108" fontSize="9" fill={residueColor(ch)} fontFamily="JetBrains Mono, monospace">
                {ch}
              </text>
            ))
          : null}
      </svg>
      <p className="hs-topbar-meta">
        Selected start {formatScientificValue(selected.start_0based ?? selected.start)} · PAM/PFS{" "}
        {formatScientificValue(selected.pam_sequence ?? selected.pfs)} · strand {formatScientificValue(selected.strand)}
      </p>
    </section>
  );
}
