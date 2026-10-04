"use client";

import { useId, useState, type ReactNode } from "react";
import { ScientificStatus } from "@/components/science/SciencePrimitives";
import { usePresentation } from "@/lib/visual/presentation";

export type FigureMeta = {
  title: string;
  shows: string;
  doesNotShow: string;
  method: string;
  xAxis?: string;
  yAxis?: string;
  units: string;
  parameters?: string;
  status?: unknown;
  source: string;
  exportName?: string;
};

export function ScientificFigure({
  meta,
  testId,
  children,
  onExportData,
}: {
  meta: FigureMeta;
  testId?: string;
  children: ReactNode;
  onExportData?: () => void;
}) {
  const [expanded, setExpanded] = useState(false);
  const dialogId = useId();
  const mode = usePresentation();
  return (
    <>
      <figure className="hs-figure" data-testid={testId}>
        <div className="hs-figure-toolbar">
          <strong>{meta.title}</strong>
          <div className="hs-actions">
            {onExportData ? (
              <button type="button" className="btn-ghost" onClick={onExportData}>
                Export data
              </button>
            ) : null}
            <button type="button" className="btn-ghost" onClick={() => setExpanded(true)} disabled={expanded}>
              Expand
            </button>
          </div>
        </div>
        {!expanded ? <div className="hs-figure-body">{children}</div> : null}
        <figcaption className="hs-figure-caption">
          <ScientificStatus status={meta.status} />
          <span>Shows: {meta.shows}</span>
          {mode === "expert" ? (
            <>
              <span>Does not show: {meta.doesNotShow}</span>
              <span>Method: {meta.method}</span>
              {meta.xAxis ? <span>X: {meta.xAxis}</span> : null}
              {meta.yAxis ? <span>Y: {meta.yAxis}</span> : null}
              <span>Units: {meta.units}</span>
              {meta.parameters ? <span>Parameters: {meta.parameters}</span> : null}
              <span>Source: {meta.source}</span>
            </>
          ) : (
            <span>Source: {meta.source}</span>
          )}
        </figcaption>
      </figure>
      {expanded ? (
        <div className="hs-figure-overlay" role="presentation" onMouseDown={() => setExpanded(false)}>
          <div
            className="hs-figure-dialog glass-medium"
            role="dialog"
            aria-modal="true"
            aria-labelledby={dialogId}
            onMouseDown={(event) => event.stopPropagation()}
          >
            <div className="hs-figure-toolbar">
              <h2 id={dialogId}>{meta.title}</h2>
              <button type="button" className="btn-secondary" onClick={() => setExpanded(false)}>
                Close
              </button>
            </div>
            <p className="hs-lede">
              {meta.shows} Status is preserved. Modebar image export (PNG/SVG) remains on the chart.
            </p>
            <div className="hs-figure-body hs-figure-body-expanded">{children}</div>
            <dl className="hs-figure-meta">
              <div>
                <dt>What it shows</dt>
                <dd>{meta.shows}</dd>
              </div>
              <div>
                <dt>What it does not show</dt>
                <dd>{meta.doesNotShow}</dd>
              </div>
              <div>
                <dt>Method</dt>
                <dd>{meta.method}</dd>
              </div>
              <div>
                <dt>Units</dt>
                <dd>{meta.units}</dd>
              </div>
              <div>
                <dt>Source</dt>
                <dd>{meta.source}</dd>
              </div>
              {meta.parameters ? (
                <div>
                  <dt>Parameters</dt>
                  <dd>{meta.parameters}</dd>
                </div>
              ) : null}
            </dl>
          </div>
        </div>
      ) : null}
    </>
  );
}
