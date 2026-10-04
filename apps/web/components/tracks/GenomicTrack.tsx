"use client";

import { useMemo, useState } from "react";
import { asList, asRecord, chartNumeric, formatScientificValue } from "@/lib/api/numeric";
import { residueColor } from "@/lib/visual/nucleotides";

type TrackFeature = {
  start: number;
  end: number;
  label: string;
  kind: string;
};

function featuresFrom(tracks: unknown, kind: string): TrackFeature[] {
  const out: TrackFeature[] = [];
  asList(tracks).forEach((track) => {
    const rec = asRecord(track);
    const name = String(rec.name || kind);
    asList(rec.features).forEach((feature) => {
      const item = asRecord(feature);
      const start = chartNumeric(item.start);
      const end = chartNumeric(item.end);
      if (start === null || end === null) return;
      out.push({
        start,
        end: end === start ? end + 1 : end,
        label: String(item.label || name),
        kind: name.toLowerCase(),
      });
    });
  });
  return out;
}

export function GenomicTrack({
  length,
  sequence,
  tracks,
  profiles,
  selected,
  onSelect,
}: {
  length: number;
  sequence?: string;
  tracks: unknown;
  profiles?: Array<Record<string, unknown>>;
  selected?: { start?: number; end?: number; kind?: string } | null;
  onSelect?: (feature: TrackFeature) => void;
}) {
  const [view, setView] = useState({ start: 0, end: 1, length: 0 });
  const features = useMemo(() => featuresFrom(tracks, "feature"), [tracks]);
  const width = 920;
  const height = 168;
  const resolved = view.length === length ? view : { start: 0, end: Math.max(length, 1), length };
  const span = Math.max(1, resolved.end - resolved.start);
  const x = (pos: number) => ((pos - resolved.start) / span) * (width - 24) + 12;

  const gc = profiles?.map((row) => ({
    mid: chartNumeric(row.midpoint ?? row.center ?? row.position),
    gc: chartNumeric(row.gc_percent ?? row.gc ?? row.gc_content),
  }));

  return (
    <div className="genomic-track" data-testid="genomic-track">
      <div className="hs-actions">
        <button type="button" className="btn-ghost" onClick={() => setView({ start: 0, end: Math.max(length, 1), length })}>
          Reset
        </button>
      </div>
      <svg
        viewBox={`0 0 ${width} ${height}`}
        role="img"
        aria-label="Linked sequence track from API coordinates"
        onWheel={(event) => {
          event.preventDefault();
          const factor = event.deltaY > 0 ? 1.15 : 1 / 1.15;
          const mid = (resolved.start + resolved.end) / 2;
          const nextSpan = Math.min(length, Math.max(40, span * factor));
          const start = Math.max(0, mid - nextSpan / 2);
          setView({ start, end: Math.min(length, start + nextSpan), length });
        }}
      >
        <rect width={width} height={height} fill="#01040B" />
        <line x1="12" x2={width - 12} y1="28" y2="28" stroke="rgba(175,203,241,0.35)" />
        {Array.from({ length: 6 }, (_, i) => {
          const pos = resolved.start + (span * i) / 5;
          return (
            <text key={i} x={x(pos)} y="18" fill="#8aabd6" fontSize="10">
              {Math.round(pos)}
            </text>
          );
        })}
        {gc?.map((row, index) => {
          if (row.mid === null || row.gc === null) return null;
          const y = 52 - (row.gc / 100) * 28;
          return <circle key={index} cx={x(row.mid)} cy={y} r="1.4" fill="#6893D0" />;
        })}
        {features.map((feature, index) => {
          const active =
            selected &&
            selected.start === feature.start &&
            selected.end === feature.end &&
            (!selected.kind || selected.kind === feature.kind);
          return (
            <rect
              key={`${feature.kind}-${index}`}
              x={x(feature.start)}
              y={feature.kind.includes("orf") ? 78 : feature.kind.includes("cpg") ? 98 : 118}
              width={Math.max(2, x(feature.end) - x(feature.start))}
              height="14"
              rx="3"
              fill={active ? "#D9E4FC" : feature.kind.includes("orf") ? "#3dd68c" : feature.kind.includes("cpg") ? "#f0c14a" : "#5b9dff"}
              tabIndex={0}
              role="button"
              aria-label={feature.label}
              onClick={() => onSelect?.(feature)}
            />
          );
        })}
        {sequence && span <= 160
          ? sequence
              .slice(Math.floor(resolved.start), Math.ceil(resolved.end))
              .split("")
              .map((ch, offset) => (
                <text
                  key={offset}
                  x={x(resolved.start + offset)}
                  y="156"
                  fontSize="9"
                  fontFamily="JetBrains Mono, monospace"
                  fill={residueColor(ch)}
                >
                  {ch}
                </text>
              ))
          : null}
      </svg>
      <p className="hs-topbar-meta">
        Coordinates are 0-based API intervals. Length {formatScientificValue(length)}. Zoom is presentation only.
      </p>
    </div>
  );
}
