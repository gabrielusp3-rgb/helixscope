"use client";

import type { ReactNode } from "react";
import { ScientificFigure } from "@/components/science/ScientificFigure";

export function InfographicFrame({
  title,
  dataSource,
  method,
  units,
  status,
  shows,
  doesNotShow,
  testId,
  children,
}: {
  title: string;
  dataSource: string;
  method: string;
  units: string;
  status?: unknown;
  shows: string;
  doesNotShow: string;
  testId?: string;
  children: ReactNode;
}) {
  return (
    <ScientificFigure
      testId={testId}
      meta={{
        title,
        shows,
        doesNotShow,
        method,
        units,
        status,
        source: dataSource,
        exportName: title,
      }}
    >
      {children}
    </ScientificFigure>
  );
}
