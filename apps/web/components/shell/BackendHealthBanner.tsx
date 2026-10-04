"use client";

import { useEffect, useState } from "react";
import { helixApi } from "@/lib/api/client";

export function BackendHealthBanner() {
  const [reachable, setReachable] = useState<boolean | null>(null);

  useEffect(() => {
    let cancelled = false;
    const ping = () => {
      helixApi
        .healthReady()
        .then(() => {
          if (!cancelled) setReachable(true);
        })
        .catch(() => {
          if (!cancelled) setReachable(false);
        });
    };
    ping();
    const timer = window.setInterval(ping, 12_000);
    return () => {
      cancelled = true;
      window.clearInterval(timer);
    };
  }, []);

  if (reachable !== false) return null;
  return (
    <div className="error-box" role="alert" data-testid="backend-unavailable">
      <strong>Backend unavailable</strong>
      <div>
        FastAPI at the configured HelixScope API origin is not reachable. This interface does not compute
        biology while the backend is down.
      </div>
    </div>
  );
}
