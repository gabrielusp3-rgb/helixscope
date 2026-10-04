"use client";

import { useEffect, useState } from "react";
import { Sidebar } from "@/components/shell/Sidebar";
import { Topbar } from "@/components/shell/Topbar";
import { CommandPalette } from "@/components/shell/CommandPalette";
import { PointerGlow } from "@/components/shell/PointerGlow";
import { BackendHealthBanner } from "@/components/shell/BackendHealthBanner";
import { ApiContractProvider } from "@/lib/api/contractIdentity";
import { ApiContractBanner } from "@/components/shell/ApiContractBanner";

export function WorkstationShell({ children }: { children: React.ReactNode }) {
  const [cmd, setCmd] = useState(false);

  useEffect(() => {
    const onKey = (event: KeyboardEvent) => {
      if ((event.ctrlKey || event.metaKey) && (event.key.toLowerCase() === "k" || event.code === "KeyK")) {
        event.preventDefault();
        event.stopPropagation();
        setCmd((open) => !open);
      }
    };
    window.addEventListener("keydown", onKey, true);
    return () => window.removeEventListener("keydown", onKey, true);
  }, []);

  return (
    <ApiContractProvider>
      <div className="hs-shell">
        <PointerGlow />
        <Sidebar />
        <div className="hs-main">
          <Topbar onCommand={() => setCmd(true)} />
          <BackendHealthBanner />
          <ApiContractBanner />
          <div className="hs-workspace" data-testid="hs-workspace">
            {children}
          </div>
        </div>
        <CommandPalette open={cmd} onClose={() => setCmd(false)} />
      </div>
    </ApiContractProvider>
  );
}
