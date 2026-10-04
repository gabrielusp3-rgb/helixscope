"use client";

import { useEffect, useRef, useState } from "react";
import { helixApi } from "@/lib/api/client";
import { asRecord } from "@/lib/api/numeric";
import { BUNDLED_MOLSTAR_IDS } from "@/components/structure/structureScene";

type PluginHandle = {
  dispose: () => void;
  managers?: { camera?: { reset?: () => void } };
  canvas3d?: { setProps?: (props: unknown) => void };
  helpers?: { viewportScreenshot?: { getImage?: (transparent?: boolean) => Promise<{ dataUri?: string }> } };
};

export function MolstarViewport({
  structureId,
  onReady,
}: {
  structureId: string;
  onReady?: (plugin: PluginHandle | null) => void;
}) {
  const host = useRef<HTMLDivElement>(null);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState("Loading molecular viewer…");

  useEffect(() => {
    const node = host.current;
    const id = structureId.toUpperCase();
    if (!BUNDLED_MOLSTAR_IDS.has(id) || !node) {
      onReady?.(null);
      return;
    }
    let disposed = false;
    let plugin: PluginHandle | null = null;
    void (async () => {
      try {
        setLoading("Fetching bundled mmCIF…");
        const { data } = await helixApi.structureCoordinateFile({ structure_id: id as "1CRN" | "1BNA" | "1RNA" });
        const mmcif = String(asRecord(data.result).mmcif || "");
        if (!mmcif) throw new Error("Bundled mmCIF was empty.");
        if (disposed) return;
        setLoading("Parsing mmCIF…");
        const { createPluginUI } = await import("molstar/lib/mol-plugin-ui");
        const { renderReact18 } = await import("molstar/lib/mol-plugin-ui/react18");
        const { DefaultPluginUISpec } = await import("molstar/lib/mol-plugin-ui/spec");
        const { Color } = await import("molstar/lib/mol-util/color");
        await import("molstar/build/viewer/theme/dark.css");
        if (disposed) return;
        const spec = {
          ...DefaultPluginUISpec(),
          layout: { initial: { isExpanded: false, showControls: false } },
          components: {
            ...(DefaultPluginUISpec().components || {}),
            controls: { left: "none", right: "none", top: "none", bottom: "none" },
            remoteState: "none",
          },
        };
        plugin = (await createPluginUI({
          target: node,
          render: renderReact18,
          spec: spec as never,
        })) as unknown as PluginHandle;
        plugin.canvas3d?.setProps?.({ renderer: { backgroundColor: Color(0x01040b) } });
        const builders = (
          plugin as unknown as {
            builders: {
              data: { rawData: (args: { data: string; label: string }) => Promise<unknown> };
              structure: {
                parseTrajectory: (data: unknown, format: string) => Promise<unknown>;
                hierarchy: { applyPreset: (traj: unknown, preset: string) => Promise<unknown> };
              };
            };
          }
        ).builders;
        setLoading("Building representation…");
        const raw = await builders.data.rawData({ data: mmcif, label: id });
        const trajectory = await builders.structure.parseTrajectory(raw, "mmcif");
        await builders.structure.hierarchy.applyPreset(trajectory, "default");
        if (disposed) {
          plugin.dispose();
          return;
        }
        setLoading("");
        onReady?.(plugin);
      } catch (err) {
        if (!disposed) setError(err instanceof Error ? err.message : "Mol* failed to initialize.");
        onReady?.(null);
      }
    })();
    return () => {
      disposed = true;
      plugin?.dispose();
      node.innerHTML = "";
    };
  }, [structureId, onReady]);

  return (
    <div className="molstar-host" data-testid="molstar-viewport">
      {loading ? <p className="hs-topbar-meta">{loading}</p> : null}
      {error ? (
        <p className="error-box" data-testid="molstar-error">
          {error} Plotly atom scatter remains available as the legacy renderer.
        </p>
      ) : null}
      <div ref={host} className="plot-host" />
    </div>
  );
}
