"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import { useRouter } from "next/navigation";
import { NAV_GROUPS, NAV_ITEMS } from "@/lib/nav";

export function CommandPalette({
  open,
  onClose,
}: {
  open: boolean;
  onClose: () => void;
}) {
  if (!open) return null;
  return <CommandPaletteOpen onClose={onClose} />;
}

function CommandPaletteOpen({ onClose }: { onClose: () => void }) {
  const router = useRouter();
  const inputRef = useRef<HTMLInputElement>(null);
  const listRef = useRef<HTMLUListElement>(null);
  const [query, setQuery] = useState("");
  const [active, setActive] = useState(0);
  const filtered = useMemo(() => {
    const q = query.trim().toLowerCase();
    if (!q) return NAV_ITEMS;
    return NAV_ITEMS.filter((item) =>
      [item.label, item.group, item.href, ...item.keywords].join(" ").toLowerCase().includes(q),
    );
  }, [query]);

  useEffect(() => {
    inputRef.current?.focus();
  }, []);

  useEffect(() => {
    const node = listRef.current?.querySelector('[data-active="true"]');
    node?.scrollIntoView({ block: "nearest" });
  }, [active, filtered]);

  useEffect(() => {
    const onKey = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        event.preventDefault();
        onClose();
      } else if (event.key === "ArrowDown") {
        event.preventDefault();
        setActive((i) => Math.min(filtered.length - 1, i + 1));
      } else if (event.key === "ArrowUp") {
        event.preventDefault();
        setActive((i) => Math.max(0, i - 1));
      } else if (event.key === "Enter") {
        event.preventDefault();
        const item = filtered[active];
        if (item) {
          router.push(item.href);
          onClose();
        }
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [filtered, active, onClose, router]);

  return (
    <div className="cmd-overlay" role="presentation" onMouseDown={onClose}>
      <div
        className="cmd-dialog glass-medium"
        role="dialog"
        aria-modal="true"
        aria-label="Module command palette"
        onMouseDown={(event) => event.stopPropagation()}
      >
        <p className="kicker">Module navigation</p>
        <input
          ref={inputRef}
          value={query}
          onChange={(event) => {
            setQuery(event.target.value);
            setActive(0);
          }}
          placeholder="Go to a workstation module"
          aria-label="Search modules"
          data-testid="command-input"
        />
        <p className="hs-topbar-meta">This palette searches HelixScope modules. It does not search scientific databases.</p>
        <ul className="cmd-list" role="listbox" ref={listRef}>
          {NAV_GROUPS.map((group) => {
            const rows = filtered.filter((item) => item.group === group);
            if (rows.length === 0) return null;
            return (
              <li key={group} className="cmd-group" role="none">
                <p className="cmd-group-label">{group}</p>
                {rows.map((item) => {
                  const index = filtered.indexOf(item);
                  return (
                    <button
                      key={item.href}
                      type="button"
                      role="option"
                      aria-selected={index === active}
                      data-active={index === active ? "true" : "false"}
                      onMouseEnter={() => setActive(index)}
                      onClick={() => {
                        router.push(item.href);
                        onClose();
                      }}
                    >
                      {item.label}
                    </button>
                  );
                })}
              </li>
            );
          })}
        </ul>
      </div>
    </div>
  );
}
