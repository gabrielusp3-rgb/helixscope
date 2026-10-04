"use client";

import { usePathname } from "next/navigation";
import { titleForPath } from "@/lib/nav";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";

export function Topbar({ onCommand }: { onCommand: () => void }) {
  const pathname = usePathname();
  const contract = useApiContract();
  const { presentationMode, setPresentationMode } = useWorkspace();
  const product = contract.productVersion;
  return (
    <header className="hs-topbar glass-thin" data-testid="hs-topbar" data-contract-status={contract.status}>
      <div>
        <p className="hs-topbar-title">{titleForPath(pathname)}</p>
        <div className="hs-topbar-meta">
          HelixScope API {contract.apiVersion || "v1"}
          {product ? ` · product ${product}` : ""}
          {contract.coreVersion ? ` · core ${contract.coreVersion}` : ""}
        </div>
      </div>
      <div className="hs-topbar-actions">
        <label className="hs-presentation-toggle">
          <span className="hs-topbar-meta">View</span>
          <select
            data-testid="presentation-mode"
            value={presentationMode}
            onChange={(event) => setPresentationMode(event.target.value as "beginner" | "expert")}
            aria-label="Beginner or expert presentation of the same result"
          >
            <option value="beginner">Beginner</option>
            <option value="expert">Expert</option>
          </select>
        </label>
        <button type="button" className="hs-cmd-btn" onClick={onCommand} data-testid="command-open">
          Pesquisar / Command
          <span className="hs-topbar-meta"> Ctrl+K</span>
        </button>
      </div>
    </header>
  );
}
