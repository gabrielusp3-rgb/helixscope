"use client";

import { useApiContract } from "@/lib/api/contractIdentity";

export function ApiContractBanner() {
  const contract = useApiContract();
  if (contract.status === "mismatch") {
    return (
      <div className="error-box" role="alert" data-testid="api-contract-mismatch">
        <strong>API/Frontend version mismatch.</strong>
        <div>
          Restart HelixScope using the supported launcher. This is a software compatibility
          state, not a scientific result. Analyses are blocked until FastAPI and this UI share
          the same OpenAPI contract.
        </div>
        <div className="hs-topbar-meta">
          expected {contract.expectedHash.slice(0, 12)} · running {contract.runningHash.slice(0, 12) || "unknown"}
        </div>
      </div>
    );
  }
  if (contract.status === "unreachable") {
    return (
      <div className="error-box" role="alert" data-testid="api-contract-unreachable">
        <strong>API identity unreachable.</strong>
        <div>
          Restart HelixScope using the supported launcher. The UI will not run scientific
          analyses until FastAPI identity is readable and matches this build.
        </div>
      </div>
    );
  }
  return null;
}
