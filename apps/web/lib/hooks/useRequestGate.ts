"use client";

import { useRef } from "react";

export function useRequestGate() {
  const generation = useRef(0);
  const controller = useRef<AbortController | null>(null);

  function next(): { signal: AbortSignal; isCurrent: () => boolean } {
    controller.current?.abort();
    controller.current = new AbortController();
    const id = ++generation.current;
    return {
      signal: controller.current.signal,
      isCurrent: () => generation.current === id,
    };
  }

  return { next };
}
