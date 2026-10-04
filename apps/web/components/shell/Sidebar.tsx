"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";
import { NAV_GROUPS, NAV_ITEMS } from "@/lib/nav";

export function Sidebar() {
  const pathname = usePathname();
  return (
    <aside className="hs-sidebar glass-thick" aria-label="Workstation navigation">
      <p className="hs-brand">HelixScope</p>
      <p className="hs-brand-sub">Scientific workstation</p>
      {NAV_GROUPS.map((group) => (
        <nav key={group} className="hs-nav-group" aria-label={group}>
          <h2>{group}</h2>
          {NAV_ITEMS.filter((item) => item.group === group).map((item) => (
            <Link
              key={item.href}
              href={item.href}
              className="hs-nav-link"
              data-active={pathname === item.href ? "true" : "false"}
            >
              {item.label}
            </Link>
          ))}
        </nav>
      ))}
    </aside>
  );
}
