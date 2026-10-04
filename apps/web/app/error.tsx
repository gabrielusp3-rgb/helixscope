"use client";

export default function RootError({
  error,
  reset,
}: {
  error: Error & { digest?: string };
  reset: () => void;
}) {
  return (
    <main className="hs-page" style={{ padding: 48 }}>
      <h1>Workstation error</h1>
      <p className="error-box">{error.message}</p>
      <button type="button" className="btn-primary" onClick={reset}>
        Retry
      </button>
    </main>
  );
}
