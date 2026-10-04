"use client";

import type { ReactNode } from "react";
import { useWorkspace, type PresentationMode } from "@/lib/workspace/WorkstationProvider";

export type { PresentationMode };

export function usePresentation(): PresentationMode {
  return useWorkspace().presentationMode;
}

export function ExpertOnly({ children }: { children: ReactNode }) {
  if (usePresentation() !== "expert") return null;
  return <>{children}</>;
}

export function BeginnerOnly({ children }: { children: ReactNode }) {
  if (usePresentation() !== "beginner") return null;
  return <>{children}</>;
}
