"use client";

import type { ReactNode } from "react";

export type AnalysisNavItem = {
  id: string;
  label: string;
  available?: boolean;
};

export function AnalysisNav({ items }: { items: AnalysisNavItem[] }) {
  const visible = items.filter((item) => item.available !== false);
  if (visible.length === 0) return null;
  return (
    <nav className="hs-analysis-nav" aria-label="Result sections" data-testid="analysis-nav">
      {visible.map((item) => (
        <a key={item.id} href={`#${item.id}`} className="hs-analysis-nav-link">
          {item.label}
        </a>
      ))}
    </nav>
  );
}

export function ResultAnchor({
  id,
  children,
}: {
  id: string;
  children: ReactNode;
}) {
  return (
    <div id={id} className="hs-result-anchor">
      {children}
    </div>
  );
}

export function AnonymousSequenceNotice() {
  return (
    <p className="hs-lede" data-testid="anonymous-sequence-notice">
      This result is sequence-derived. Organism, gene, chromosome and function are unknown
      until a separate NCBI or BLAST retrieval provides evidence. HelixScope does not guess
      identity from an anonymous paste.
    </p>
  );
}

export function RemotePrivacyNotice({ service, what }: { service: string; what: string }) {
  return (
    <p className="hs-lede" data-testid="remote-privacy-notice">
      Remote transmission: {what} may be sent to {service}. This is not a silent background
      upload. Secrets and local file paths stay on this workstation.
    </p>
  );
}
