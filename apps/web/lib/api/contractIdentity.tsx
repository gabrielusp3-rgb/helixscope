"use client";

import { createContext, useContext, useEffect, useMemo, useState, type ReactNode } from "react";
import { helixApi } from "@/lib/api/client";
import { asRecord } from "@/lib/api/numeric";
import { EXPECTED_OPENAPI_SHA256 } from "@/lib/api/generated/contract";
import { evaluateApiContract } from "@/lib/api/contractEvaluate";

export type ContractStatus = "pending" | "ok" | "mismatch" | "unreachable";

export type ContractValue = {
  status: ContractStatus;
  scienceLocked: boolean;
  productVersion: string;
  coreVersion: string;
  apiVersion: string;
  expectedHash: string;
  runningHash: string;
};

export { evaluateApiContract } from "@/lib/api/contractEvaluate";

const EMPTY: ContractValue = {
  status: "pending",
  scienceLocked: true,
  productVersion: "",
  coreVersion: "",
  apiVersion: "",
  expectedHash: EXPECTED_OPENAPI_SHA256,
  runningHash: "",
};

const ContractContext = createContext<ContractValue>(EMPTY);

export function ApiContractProvider({ children }: { children: ReactNode }) {
  const [state, setState] = useState<ContractValue>(EMPTY);

  useEffect(() => {
    let cancelled = false;
    const controller = new AbortController();
    const timer = window.setTimeout(() => controller.abort(), 8000);
    helixApi
      .systemIdentity(controller.signal)
      .then(({ data }) => {
        if (cancelled) return;
        const rec = asRecord(data);
        const runningHash = String(rec.openapi_sha256 || "");
        const verdict = evaluateApiContract(EXPECTED_OPENAPI_SHA256, runningHash);
        setState({
          status: verdict.status,
          scienceLocked: verdict.scienceLocked,
          productVersion: String(rec.product_version || ""),
          coreVersion: String(rec.core_version || ""),
          apiVersion: String(rec.api_version || ""),
          expectedHash: EXPECTED_OPENAPI_SHA256,
          runningHash,
        });
      })
      .catch(() => {
        if (!cancelled) {
          setState({
            ...EMPTY,
            status: "unreachable",
            scienceLocked: true,
          });
        }
      })
      .finally(() => {
        window.clearTimeout(timer);
      });
    return () => {
      cancelled = true;
      controller.abort();
      window.clearTimeout(timer);
    };
  }, []);

  const value = useMemo(() => state, [state]);
  return <ContractContext.Provider value={value}>{children}</ContractContext.Provider>;
}

export function useApiContract(): ContractValue {
  return useContext(ContractContext);
}
