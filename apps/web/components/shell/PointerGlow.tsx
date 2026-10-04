"use client";

import { useEffect } from "react";

export function PointerGlow() {
  useEffect(() => {
    const onMove = (event: PointerEvent) => {
      const root = document.documentElement;
      root.style.setProperty("--mx", `${(event.clientX / window.innerWidth) * 100}%`);
      root.style.setProperty("--my", `${(event.clientY / window.innerHeight) * 100}%`);
    };
    window.addEventListener("pointermove", onMove, { passive: true });
    return () => window.removeEventListener("pointermove", onMove);
  }, []);
  return null;
}
