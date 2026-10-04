import Link from "next/link";

export default function NotFound() {
  return (
    <main className="hs-page" style={{ padding: 48 }}>
      <p className="kicker">HelixScope</p>
      <h1>Module not found</h1>
      <p className="hs-lede">This route is not a workstation module.</p>
      <Link href="/overview" className="btn-primary" style={{ display: "inline-block" }}>
        Overview
      </Link>
    </main>
  );
}
