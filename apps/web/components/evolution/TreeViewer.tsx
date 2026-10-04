"use client";

import { useMemo, useState } from "react";
import { asList, asRecord, formatScientificValue } from "@/lib/api/numeric";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";

const TAX_COLORS = ["#3dd68c", "#5b9dff", "#f0c14a", "#ff6b8a", "#c084fc", "#5eead4", "#91B4E4"];

export function TreeViewer({ result }: { result: Record<string, unknown> }) {
  const layouts = asRecord(result.layouts);
  const rectangular = asRecord(layouts.rectangular);
  const nodes = asList(rectangular.nodes).map((n) => asRecord(n));
  const edges = asList(rectangular.edges).map((e) => asRecord(e));
  const leavesMeta = asList(result.leaves).map((item) => asRecord(item));
  const taxonomyLayer = asRecord(result.taxonomy_layer);
  const [zoom, setZoom] = useState(1);
  const [pan, setPan] = useState({ x: 40, y: 24 });
  const [query, setQuery] = useState("");
  const { setSelectedTreeLeaf, selectedTreeLeaf, setSelectedMsaId } = useWorkspace();

  const taxColor = useMemo(() => {
    const map = new Map<string, string>();
    const groups = new Map<string, string>();
    let i = 0;
    for (const leaf of leavesMeta) {
      const tax = asRecord(leaf.taxonomy);
      const group = String(tax.rank || tax.scientific_name || leaf.organism_display || "unannotated");
      if (!groups.has(group)) {
        groups.set(group, TAX_COLORS[i % TAX_COLORS.length]);
        i += 1;
      }
      const id = String(leaf.tree_id || leaf.identifier || leaf.id || "");
      if (id) map.set(id, groups.get(group) || TAX_COLORS[0]);
    }
    return map;
  }, [leavesMeta]);

  let maxX = 1;
  let maxY = 1;
  for (const node of nodes) {
    maxX = Math.max(maxX, Number(node.x) || 0);
    maxY = Math.max(maxY, Number(node.y) || 0);
  }
  const bounds = { maxX: maxX || 1, maxY: maxY || 1 };
  const width = 900;
  const height = 420;
  const sx = ((width - 160) / bounds.maxX) * zoom;
  const sy = ((height - 80) / Math.max(1, bounds.maxY)) * zoom;
  const needle = query.trim().toLowerCase();

  if (nodes.length === 0) {
    return <p>No API tree layout is present.</p>;
  }

  function selectLeaf(label: string) {
    setSelectedTreeLeaf(label);
    const hit = leavesMeta.find((item) => String(item.tree_id) === label || String(item.identifier) === label);
    if (hit && hit.index !== undefined) setSelectedMsaId(String(hit.identifier ?? label));
  }

  return (
    <div>
      <div className="hs-actions">
        <button type="button" className="btn-secondary" onClick={() => setZoom((z) => Math.min(4, z * 1.15))}>
          Zoom in
        </button>
        <button type="button" className="btn-secondary" onClick={() => setZoom((z) => Math.max(0.4, z / 1.15))}>
          Zoom out
        </button>
        <button
          type="button"
          className="btn-ghost"
          onClick={() => {
            setZoom(1);
            setPan({ x: 40, y: 24 });
          }}
        >
          Fit
        </button>
        <label>
          Search leaf
          <input
            type="search"
            value={query}
            onChange={(event) => setQuery(event.target.value)}
            data-testid="tree-leaf-search"
          />
        </label>
      </div>
      {taxonomyLayer.status ? (
        <p className="hs-topbar-meta" data-testid="taxonomy-annotation">
          TAXONOMIC ANNOTATION · {String(taxonomyLayer.status)} · {String(taxonomyLayer.note || "does not change topology")}
        </p>
      ) : null}
      <svg
        className="tree-svg"
        viewBox={`0 0 ${width} ${height}`}
        data-testid="tree-viewer"
        role="img"
        aria-label="Phylogenetic tree from API layout"
      >
        <g transform={`translate(${pan.x} ${pan.y})`}>
          {edges.map((edge, index) => (
            <path
              key={index}
              d={`M ${Number(edge.x0) * sx} ${Number(edge.y0) * sy} V ${Number(edge.y1) * sy} H ${Number(edge.x1) * sx}`}
              fill="none"
              stroke="#6893D0"
              strokeWidth="1.2"
            />
          ))}
          {nodes.map((node) => {
            const x = Number(node.x) * sx;
            const y = Number(node.y) * sy;
            const leaf = Boolean(node.is_leaf);
            const label = String(node.label);
            const selected = selectedTreeLeaf === label || (needle && label.toLowerCase().includes(needle));
            const fill = leaf
              ? taxColor.get(label) || (selected ? "#D9E4FC" : "#91B4E4")
              : "#2355A0";
            return (
              <g key={String(node.id)}>
                <circle
                  cx={x}
                  cy={y}
                  r={leaf ? 5 : 3}
                  fill={selected && leaf ? "#D9E4FC" : fill}
                  tabIndex={0}
                  role="button"
                  aria-label={label}
                  onClick={() => {
                    if (!leaf) return;
                    selectLeaf(label);
                  }}
                />
                {leaf ? (
                  <text x={x + 8} y={y + 4} fill="#D9E4FC" fontSize="12">
                    {label}
                  </text>
                ) : null}
              </g>
            );
          })}
        </g>
      </svg>
      <p className="hs-topbar-meta">
        Layout coordinates are from the API. Pan/zoom is camera-like presentation only. Support{" "}
        {formatScientificValue(asRecord(result.support).status ?? result.support)}. Absence of support is
        N/A, not 0. Leaf color from taxonomy is TAXONOMIC ANNOTATION, not phylogenetic inference.
      </p>
    </div>
  );
}
